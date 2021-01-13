This folder includes several lists for hooks.

- Commit hook:

    - `absolute_path_list`: CBIG repo avoids using absolute paths. The list includes a set of keywords that might appear in invalid absolute paths. If any of these keywords are detected, the commit will be aborted. 

    - `banned_keyword_list`: Currently, `cp` is not allowed in CBIG repo. The list includes banned keywords. If any of these keywords are detected, the commit will be aborted. 

- Push hook:

    - `exclude_list`: All projects in this list will not be checked.

    - `exclude_replication_list`: For projects in this list, we will skip the part of checking replication folder.

    - `exclude_example_list`: For projects in this list, we will skip the part of checking example folder.

    - `exclude_standalone_list`: For projects in this list, we will skip the part of checking if the project has a standalone script.
