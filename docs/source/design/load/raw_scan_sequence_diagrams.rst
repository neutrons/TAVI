Raw Scan Load Sequence Diagram
==============================

Main Flow
---------

```mermaid
sequenceDiagram
    participant EventBroker
    participant ProjectModel
    participant RawScanLoadController
    participant FileStore
    participant RawScanClassifier
    participant LoaderRegistry
    participant ORNLSpiceLoader
    LoaderRegistry ->> LoaderRegistry: init(FileStore)
    ORNLSpiceLoader -->> LoaderRegistry: is registered via main
    LoaderRegistry ->> LoaderRegistry: init loader instance with appropriate FileStore
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
        RawScanClassifier -->> RawScanLoadController: (ORNL SPICE) stubbed

        RawScanLoadController ->> LoaderRegistry: buildLoader(classification, DataService)
        LoaderRegistry -->> RawScanLoadController: ORNLSpiceLoader
        RawScanLoadController ->> ORNLSpiceLoader: load(raw scan file (handler?))

        ORNLSpiceLoader ->> ORNLSpiceLoader: setup
        ORNLSpiceLoader ->> ORNLSpiceLoader: parse metadata
        ORNLSpiceLoader ->> ORNLSpiceLoader: parse ub
        ORNLSpiceLoader ->> ORNLSpiceLoader: parse scan values
        ORNLSpiceLoader ->> FileStore: fetch external metadata files
        FileStore -->> ORNLSpiceLoader: return external metadata or None
        ORNLSpiceLoader ->> ORNLSpiceLoader: parse external metadata
        ORNLSpiceLoader ->> ORNLSpiceLoader: adapt parsed data into RawScan
        ORNLSpiceLoader -->> RawScanLoadController: RawScan
        RawScanLoadController ->>RawScanLoadController: append to result List
    end
    RawScanLoadController -->>ProjectModel: List[RawScan]
    ProjectModel ->> ProjectModel : update TaviData.scans
    ProjectModel ->> EventBroker : emit RawScanListUpdateEvent
```

Classificatiion Flow
--------------------

```mermaid
sequenceDiagram
	participant RawScanLoadController
	participant RawScanClassifier
	participant LoaderRegistry
    participant Loader
    participant ORNLSpiceLoader
    participant RuleBasedClassifier
	participant RuleSet
	participant Rule
    ORNLSpiceLoader ->> LoaderRegistry: Register with Registry via main
    LoaderRegistry ->> LoaderRegistry: Init loaders with appropriate FileStore
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
		RawScanClassifier ->> RawScanClassifier:return best match
		RawScanClassifier -->> RawScanLoadController: (ORNL SPICE) stubbed
		RawScanLoadController ->>RawScanLoadController: ...
		RawScanLoadController ->>RawScanLoadController: append to result List
	RawScanLoadController ->>RawScanLoadController: append to result List
	end
```

Disk Access Flow
----------------

```mermaid
sequenceDiagram
	participant RawScanLoadController
	participant DiskServiceInterface
	participant LocalDiskService
	RawScanLoadController ->>LocalDiskService: get scan(s) at (path)
	LocalDiskService ->>pathlib: get list of ALL files @ path from disk
	loop foreach filepath
		LocalDiskService ->> LocalDiskService: get file size(path)
		LocalDiskService ->> pathlib:.stat()
		pathlib -->> LocalDiskService : filesize
			alt 0 < filesize < threshold
				LocalDiskService ->>LocalDiskService: skip
			else
				LocalDiskService ->>LocalDiskService: append to result list
			end
	end
	LocalDiskService -->> RawScanLoadController: List[str]
```
