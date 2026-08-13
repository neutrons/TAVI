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
        ProjectModel ->> ProjectModel : update TaviData.raw_scans
        loop foreach RawScan
            ProjectModel ->> EventBroker : publish RawScanAppendEvent
        end


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
