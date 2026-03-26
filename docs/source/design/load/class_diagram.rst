Class Diagram
=============
.. mermaid::

    classDiagram
        class UUID {
            MD5Sum value
        }
        RawScanLoadControllerInterface : +load_raw_scans(String path)
        class RawScanLoadController {
            +FileStore filestore
            +LoaderRegistry loader_registry
            +RawScanClassifier raw_scan_classifier
            +load_scan(String path, Loader loader=None)
            +load_folder(String path, Loader loader=None, quick=Config[backend.load.quick])
            +load_scans(List[String] paths, Loader loader=None, quick=Config[backend.load.quick])
        }
        class FileStoreInterface {
            +fetch_files_at(String path)
            -validate_file(Path path)
        }
        class LocalFileStore {
            +fetch_files_at(String path)
            -validate_file(Path path)
            -read_file(Path path)
            -get_file_ext(Path path)
            -get_file_size_mb(Path path)
        }
        class RawScanClassifier{
            LoaderRegistry loader_registry
            +classify(str path) : RawScanType
        }
        class RawScanType~StrEnum~{
            ORNLSpice
            UNKNOWN
        }
        class LoaderRegistry {
            +getLoaders(FileStoreInterface filestore) : List~Loader~
            +buildLoader(RawScanType type, FileStoreInterface filestore) : Loader
        }
        class Loader{
            +__init__(FileStoreInterface filestore)
            +load(str path) : Scan
            +get_score(str path): int
            +generate_uuid(str path) : UUID
        }
        class ORNLSpiceLoader{
            +get_score(str path) : int
        }
        class RuleBasedClassifier{
            RuleSet rule_set
            +get_score(str path) : int
            -get_score_from_rules(str path) : int
            -get_score_from_rule(str path, Rule rule) : int
        }
        class RuleSet{
            List~Rule~ Rules
            +get_rules():List~Rule~
            +add_rule(Rule rule)
        }
        class RuleInterface{
            +get_score(str path, FileStoreInterface filestore) : int
        }
        class Scan {
            ScanValues values
            RawScanMetadata metadata
            TaviMetadata: tavi_metadata
            Provenance provenance
            UUID uuid
        }
        class RawScan {
        }
        class Provenance {
            Str: raw_file_path
            Dict~UUID,int~ contributing_scans
        }
        class RawScanMetadata{
            Dict~str, Any~ data
            LoaderEnum loader_name
        }
        class TaviMetadata{
            Tuple~str, str~ default_axis_cols
            str normalization
        }
        class ScanValues {
            Dict~str, list~ data
        }
        class Fit
        class ProjectModel {
            TaviData project
            +load_raw_scans(List~String~ paths)
            +load_raw_scan_folder(String path)
        }
        class TaviData {
            Dict~UUID,RawScan~ raw_scans
            Dict~UUID,ComboScan~ combo_scans
            Dict~UUID,Fit~ fits
            Dict~UUID,Plot~ plots
        }

        RawScanLoadControllerInterface <|-- RawScanLoadController
        RawScanLoadController *-- FileStoreInterface
        RawScanLoadController *-- LoaderRegistry
        RawScanLoadController *-- RawScanClassifier

        RawScanClassifier *-- LoaderRegistry
        Loader *-- FileStoreInterface
        Loader <|-- ORNLSpiceLoader
        LoaderRegistry o-- Loader

        FileStoreInterface <|-- LocalFileStore

        RuleSet o-- RuleInterface
        RuleInterface <|-- Rule
        RuleBasedClassifier *-- RuleSet
        ORNLSpiceLoader *-- RuleBasedClassifier

        Scan <|-- RawScan
        Scan <|-- ComboScan
        Scan *-- RawScanMetadata
        Scan *-- ScanValues
        Scan *-- Provenance
        Scan *-- UUID
        Scan *-- TaviMetadata
        TaviData o-- Fit
        TaviData o-- RawScan
        TaviData o-- ComboScan
        TaviData o-- Plot
        TaviData o-- UUID
        ProjectModel *-- TaviData
        ProjectModel *-- RawScanLoadControllerInterface
