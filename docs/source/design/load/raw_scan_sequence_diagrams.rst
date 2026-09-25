Raw Scan Load Sequence Diagram
==============================

Main Flow
---------

.. mermaid::

    sequenceDiagram
        participant EventBroker
        participant ProjectModel
        participant RawScanLoadController
        participant FileStore
        participant RawScanClassifier
        participant LoaderRegistry
        participant ORNLSpiceLoader
        LoaderRegistry ->> LoaderRegistry: init() registers built-in loaders
        ORNLSpiceLoader -->> LoaderRegistry: registered under its get_scan_type()
        LoaderRegistry ->> LoaderRegistry: set_filestore(FileStore) via library.init()
        ProjectModel ->>+RawScanLoadController: load raw scans (path)
        RawScanLoadController ->>FileStore: get files(s) at (path)
        FileStore ->>FileStore: validate file on disk
        FileStore -->> RawScanLoadController: Absolute File Paths(s)
        loop foreach File
            RawScanLoadController ->> RawScanClassifier: classify input schema
            RawScanClassifier ->> LoaderRegistry: getLoaders
            loop foreach Loader:
                RawScanClassifier ->> ORNLSpiceLoader: getScore(str path)
                RawScanClassifier ->> RawScanClassifier: update best score
            end
            RawScanClassifier -->> RawScanLoadController: RawScanType.ORNLSpice

            RawScanLoadController ->> LoaderRegistry: get_loader(classification)
            LoaderRegistry -->> RawScanLoadController: ORNLSpiceLoader
            RawScanLoadController ->> ORNLSpiceLoader: load(file_path)

            ORNLSpiceLoader ->> ORNLSpiceLoader: generate_uuid (md5 of file text)
            ORNLSpiceLoader ->> ORNLSpiceLoader: parse_scan_values
            ORNLSpiceLoader ->> ORNLSpiceLoader: parse_metadata
            ORNLSpiceLoader ->> ORNLSpiceLoader: parse_tavi_metadata
            ORNLSpiceLoader ->> ORNLSpiceLoader: create_provenance
            ORNLSpiceLoader ->> FileStore: read sibling UBConf file
            FileStore -->> ORNLSpiceLoader: return UB configuration text
            ORNLSpiceLoader ->> ORNLSpiceLoader: parse_external_metadata, merge into meta.data
            ORNLSpiceLoader ->> ORNLSpiceLoader: adapt_scan_data -> RawScan
            ORNLSpiceLoader -->> RawScanLoadController: RawScan
            RawScanLoadController ->>RawScanLoadController: append to result List
        end
        RawScanLoadController -->>ProjectModel: List[RawScan]
        loop foreach RawScan
            alt uuid already in TaviData.raw_scans
                ProjectModel ->> ProjectModel : skip - keep the stored scan, publish nothing
            else
                ProjectModel ->> ProjectModel : update TaviData.raw_scans
                ProjectModel ->> EventBroker : publish RawScanAppendEvent
            end
        end


Re-loading a folder
-------------------

A folder may be loaded as often as the user likes, including while an
experiment is still writing to it. ``ProjectModel.load_raw_scan_from_folder``
decides per scan whether it is new, keyed on the uuid the loader derived from
the **file's text** (``generate_uuid``, md5):

- **Unchanged file** — same text, same uuid. The scan is skipped: the copy
  already in ``TaviData.raw_scans`` is kept rather than overwritten (which
  would discard edits to its writable ``tavimeta``), and no
  ``RawScanAppendEvent`` is published.
- **File appended to since the last load** — different text, different uuid.
  It loads as a *new* scan and is announced normally, so the earlier, shorter
  scan stays in the project alongside it; plots and fits built on it keep
  resolving.

Deduplicating here, rather than in the view, is what lets
``TreeViewWidget.add_item_at_path`` keep its duplicate-uuid guard as a real
invariant: ``uuid_map`` holds one item per uuid, so a second insert of the same
uuid would leave the first row orphaned in the tree. A uuid reaching the tree
twice now means a bug upstream, not a re-load, and it is still raised as such.


Classification Flow
--------------------

.. mermaid::

    sequenceDiagram
        participant RawScanLoadController
        participant RawScanClassifier
        participant LoaderRegistry
        participant Loader
        participant ORNLSpiceLoader
        participant RuleBasedClassifier
        participant RuleSet
        participant Rule
        ORNLSpiceLoader ->> LoaderRegistry: Registered in LoaderRegistry.__init__
        LoaderRegistry ->> LoaderRegistry: set_filestore propagates FileStore to all loaders
        loop foreach File
            RawScanLoadController ->> RawScanClassifier: classify input schema
            RawScanClassifier ->> LoaderRegistry : get loaders
            loop foreach Loader
            RawScanClassifier ->> Loader: generate classification score
            Loader ->> ORNLSpiceLoader: get score
            ORNLSpiceLoader ->> RuleBasedClassifier : get score
                loop foreach Rule in RuleSet
                    RuleBasedClassifier ->> Rule: does this match?
                    Rule -->> RuleBasedClassifier: result
                    RuleBasedClassifier ->> RuleBasedClassifier: append to results List
                end
                RuleBasedClassifier ->> RuleBasedClassifier: Calculate score
                RuleBasedClassifier -->> ORNLSpiceLoader: score
                ORNLSpiceLoader -->> RawScanClassifier: score
            end
            RawScanClassifier ->> RawScanClassifier: keep highest score (ties keep the first)
            RawScanClassifier -->> RawScanLoadController: winning loader's RawScanType
            RawScanLoadController ->>RawScanLoadController: load via that loader
            RawScanLoadController ->>RawScanLoadController: append to result List
        end

Note: ``RawScanClassifier`` starts from ``(RawScanType.NONE, 0)`` and only
replaces it on a **strictly** higher score, so a file every loader rejects
classifies as ``RawScanType.NONE`` and resolves to ``DefaultLoader``, which raises
on ``load()``.


Disk Access Flow
----------------

.. mermaid::

    sequenceDiagram
        participant RawScanLoadController
        participant LocalFileStore
        participant pathlib
        RawScanLoadController ->>LocalFileStore: fetch_files_at(path)
        LocalFileStore ->>LocalFileStore: raise RuntimeError if path missing or is a file
        LocalFileStore ->>pathlib: iterdir() - all entries at path
        loop foreach entry
            LocalFileStore ->> LocalFileStore: validate_file(path)
            LocalFileStore ->> LocalFileStore: _is_real_file(path)
            LocalFileStore ->> pathlib: .stat().st_size
            pathlib -->> LocalFileStore : filesize
                alt not a real file OR size > library.filestore.raw.size-limit
                    LocalFileStore ->>LocalFileStore: skip
                else
                    LocalFileStore ->>LocalFileStore: append absolute path to result list
                end
        end
        LocalFileStore ->>LocalFileStore: sort()
        LocalFileStore -->> RawScanLoadController: List[str]

The size limit is a **maximum**, not a minimum: files *larger* than
``library.filestore.raw.size-limit`` (1 MB by default) are skipped. Directory
entries are skipped too, so ``fetch_files_at`` never recurses.
