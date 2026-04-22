Class Diagram
=============
.. mermaid::

    classDiagram
        class UUID {
            String value
        }
        RawScanLoadControllerInterface : +load_raw_scans(String path)
        class RawScanLoadController {
            +FileStore filestore
            +LoaderRegistry loader_registry
            +RawScanClassifier raw_scan_classifier
            -lookup_loader(String file_path) : AbstractLoader
            +__init__(FileStoreInterface filestore)
            +load_file(String file_path, LoaderInterface loader=None) : RawScan
            +load_folder(String folder_path, Loader loader=None, quick=Config[library.storage.raw.classification.quick]) : List~RawScan~
            +load_files(List[String] file_paths, LoaderInterface loader=None, quick=Config[library.storage.raw.classification.quick]) : List~RawScan~
        }
        class FileStoreInterface {
            +fetch_files_at(String path) : list~String~
            +validate_file(String file_path) : bool
            +write_user_data_file(String file_subpath, String value)
            +write_text_file(String file_path, String value)
            +read_text_file(String file_path) : String
            +get_file_ext(String file_path) : String
            +get_file_size_mb(String file_path) : float
        }
        class LocalFileStore {
            -is_real_file(Path file_path, bool throws=False) : bool
        }
        class RawScanClassifier{
            LoaderRegistry loader_registry
            +get_classification(str file_path) : RawScanType
        }
        class RawScanType~StrEnum~{
            ORNLSpice
            UNKNOWN
        }
        class LoaderRegistry {
            -register_loader(AbstractLoader loader)
            -refresh_filestore()
            +__init__()
            +register(String key, AbstractLoader loader)
            +set_filestore(FileStoreInterface filestore)
            +get_loaders() : List~AbstractLoader~
            +get_loader(String key) : AbstractLoader
        }
        class LoaderInterface{
            +load(String path) : Scan
            +get_scan_type() : RawScanType
            +get_score(String path): int
            +parse_metadata(String path) : ScanMetadata
            +parse_tavi_metadata(String path) : TaviMetadata
            +parse_scan_values(String path) : ScanData
            +parse_external_metadata(String path) : dict~string, Any~
            +adapt_scan_data(ScanMetadata meta, ScanData values, TaviMetadata tavi_meta) : RawScan
        }
        class AbstractLoader{
            FileStoreInterface filestore
            +__init__(FileStoreInterface filestore)
            +generate_uuid(str file_path) : UUID
        }
        class ORNLSpiceLoader{
            RuleBasedClassifier classifier
            ORNLSpiceRuleSet classification_rules
        }
        class RuleBasedClassifier{
            RuleSet rule_set
            +get_score(str path, RuleSet rule_set) : int
            -get_score_from_rules(str path, RuleSet rule_set) : int
        }
        class ORNLSpiceRuleSet{
        }
        class RuleSet{
            dict~RuleInterface, int~ rules
            +register(RuleInterface rule, int weight)
            +validate()
            +get_rules():List~RuleInterface~
            +get_weight(RuleInterface rule) : int
        }
        class RuleInterface{
            +get_score(str path, FileStoreInterface filestore) : int
        }
        class Scan {
            UUID uuid
            ScanData data
            ScanMetadata metadata
            TaviMetadata: tavimeta
            Provenance prov
        }
        class RawScan {
        }
        class Provenance {
            Str: raw_file_path
            Dict~UUID,int~ contributing_scans
        }
        class ScanMetadata{
            Dict~str, Any~ data
            LoaderEnum loader_name
        }
        class TaviMetadata{
            Tuple~str, str~ default_axis_cols
            String normalization = None
            String friendly_name
            String friendly_path
        }
        class ScanData {
            Dict~str, list~ data
        }

        RawScanLoadControllerInterface <|-- RawScanLoadController
        RawScanLoadController *-- FileStoreInterface
        RawScanLoadController *-- LoaderRegistry
        RawScanLoadController *-- RawScanClassifier

        RawScanClassifier *-- LoaderRegistry
        AbstractLoader *-- FileStoreInterface
        AbstractLoader <|-- ORNLSpiceLoader
        LoaderRegistry o-- AbstractLoader

        FileStoreInterface <|-- LocalFileStore

        RuleSet o-- RuleInterface
        RuleInterface <|-- Rule
        RuleBasedClassifier *-- RuleSet
        ORNLSpiceLoader *-- RuleBasedClassifier

        Scan <|-- RawScan
        Scan <|-- ComboScan
        Scan *-- ScanMetadata
        Scan *-- ScanData
        Scan *-- Provenance
        Scan *-- UUID
        Scan *-- TaviMetadata


.. TaviData o-- Fit
.. TaviData o-- RawScan
.. TaviData o-- ComboScan
.. TaviData o-- Plot
.. TaviData o-- UUID
.. ProjectModel *-- TaviData
.. ProjectModel *-- RawScanLoadControllerInterface

.. class ProjectModel {
..     TaviData project
..     +load_raw_scans(List~String~ paths)
..     +load_raw_scan_folder(String path)
.. }
.. class TaviData {
..     Dict~UUID,RawScan~ raw_scans
..     Dict~UUID,ComboScan~ combo_scans
..     Dict~UUID,Fit~ fits
..     Dict~UUID,Plot~ plots
.. }
