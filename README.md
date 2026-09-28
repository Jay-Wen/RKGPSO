# RKG-BOPSO 



1. Open the `bin` folder and double-click `rkgpso.exe` to run the program.

2. The program will automatically read the `run.ini` file located in the `code` directory.
3. If the program is moved or the project is changed, update the following two paths in `run.ini`:

```
[Run]
ProjectDirectory=E:\code\project_example_shifted
InitialDatabase=E:\code\project_example_shifted\DATA\example.mdb
```

- `ProjectDirectory`: The absolute path to the project folder.
- `InitialDatabase`: The absolute path to the initial alignment scheme `.mdb` file.

Use absolute paths. Do not add quotation marks around the paths on the right-hand side of the equals sign (`=`).

The calculation results are saved in the `RESULT` and `DATA` folders of the project directory.
