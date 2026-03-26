Workflow Diagram
================
.. mermaid::
    block-beta
        columns 13
        space space space space space space space space X["Parse Extra Metadata"] space space  space space
        space space space space space space space space space space space  space space
        space space space space space space space space Y["Read UB conf File"] space space space space
        space space space space space space space space space space space space space
        space space 1["File Path"] space space space 2["FileStore"] space Z["Parse Metadata/Scan Values"] space space space space
        space space space space space space space space space space space space space
        S(("start")) space A["Validate File Size"] space B["Classify File Format"] space C["Build Loader"] space G{"If ORNL Spice"} space W["Adapt to Scan Object"] space V(("end"))
        space space space space space space space space space space space space space
        space space F["File Path"] space E["File Classification"] space D["ORNLSpiceLoader"] space  space space H["Scan Object"] space space

        1 --> A
        S --> A
        A --> B
        B --> C
        A --> F
        F --> B
        B --> E
        E --> C
        2 --> C
        C --> D
        C --> G
        G -- "yes"--> Z
        G --> W
        Z --> Y
        Y --> X
        W --> V
        W --> H
        X --> W
```
