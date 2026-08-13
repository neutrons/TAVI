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
            +read_user_data_file(String file_subpath) : String
            +write_text_file(String file_path, String value)
            +read_text_file(String file_path) : String
            +get_file_ext(String file_path) : String
            +get_file_size_mb(String file_path) : float
            +get_file_name(String file_path) : String
            +get_parent(String file_path) : String
            +join_path(String root_path, String target_path) : String
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
            NONE
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
            +get_score(String path): float
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
        class DefaultLoader{
            +get_scan_type() : RawScanType.NONE
            +get_score(str path) : 0
        }
        class RuleBasedClassifier{
            FileStoreInterface filestore
            +__init__(FileStoreInterface filestore)
            +set_filestore(FileStoreInterface filestore)
            +get_score(str path, RuleSet rule_set) : float
            -get_score_from_rules(str path, RuleSet rule_set) : float
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
            Str raw_file
            Dict~UUID,int~ contributing_scans
        }
        class ScanMetadata{
            Dict~str, Any~ data
            Dict~str, List~str~~ categories
            +by_category() Dict~str, Any~
        }
        class TaviMetadata{
            Tuple~str, str~ default_axis
            String friendly_name
            String friendly_path
            Optional~Tuple~str, float~~ normalization = None
        }
        class ScanData {
            Dict~str, list~float~~ data
        }

        RawScanLoadControllerInterface <|-- RawScanLoadController
        RawScanLoadController *-- FileStoreInterface
        RawScanLoadController *-- LoaderRegistry
        RawScanLoadController *-- RawScanClassifier

        RawScanClassifier *-- LoaderRegistry
        AbstractLoader *-- FileStoreInterface
        AbstractLoader <|-- ORNLSpiceLoader
        AbstractLoader <|-- DefaultLoader
        LoaderInterface <|-- AbstractLoader
        LoaderRegistry o-- AbstractLoader

        FileStoreInterface <|-- LocalFileStore

        RuleSet <|-- ORNLSpiceRuleSet
        RuleSet o-- RuleInterface
        RuleBasedClassifier ..> RuleSet : scores against
        ORNLSpiceLoader *-- RuleBasedClassifier
        ORNLSpiceLoader *-- ORNLSpiceRuleSet

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
