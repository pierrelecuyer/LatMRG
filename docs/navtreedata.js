/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "LatMRG Guide", "index.html", [
    [ "Overview", "index.html#overview", [
      [ "Content Outline", "index.html#latmrg_guide_out", null ],
      [ "Presentation of LatMRG", "index.html#presentation", null ]
    ] ],
    [ "Theory of Random Numbers", "background.html", [
      [ "Linear Congruential Engines", "background.html#back_gen", [
        [ "Simple Congruential Generators", "background.html#back_gen_gen", null ],
        [ "Carry Generators", "background.html#back_gen_carry", null ],
        [ "Combined Generators", "background.html#back_gen_comb", null ],
        [ "Matrix Generators", "background.html#back_gen_mat", null ]
      ] ],
      [ "Lattices and Merit", "background.html#back_other", [
        [ "Lattice of a MRG", "background.html#back_lat", null ],
        [ "Measures on MRGs", "background.html#back_merit", [
          [ "Period Length", "background.html#back_merit_per", null ],
          [ "Figures of Merit", "background.html#back_merit_spectral", null ]
        ] ],
        [ "Minkowski reduced basis", "background.html#minkowski", null ]
      ] ],
      [ "A Good Parameter Choice", "background.html#back_param", null ],
      [ "Large numbers, matrices and polynomials", "background.html#numbers_sec", null ]
    ] ],
    [ "Compilation and Executables Usage", "usage.html", [
      [ "LatMRG executable", "usage.html#usage_exec", [
        [ "Searching m and k parameters", "usage.html#usage_exec_mk", null ],
        [ "Checking for full period", "usage.html#usage_exec_period", null ],
        [ "Testing lattice structure for MRG", "usage.html#usage_exec_lattest", null ],
        [ "Searching for new MRG", "usage.html#usage_exec_search", null ],
        [ "Usage of the executable", "usage.html#usage_exec_usage", null ]
      ] ],
      [ "Building LatMRG", "usage.html#compilation", null ]
    ] ],
    [ "Configuration Files Synthax and Tags", "conf_file.html", [
      [ "A short note on XML", "conf_file.html#xml_note", null ],
      [ "Program wide tags", "conf_file.html#conf_all", null ],
      [ "Search for m and k", "conf_file.html#conf_mk", null ],
      [ "Full period test", "conf_file.html#conf_period", null ],
      [ "Generator test and search", "conf_file.html#conf_lattest_seek", null ]
    ] ],
    [ "Tutorial and Library Usage", "tutorial.html", [
      [ "Classes of LatMRG", "tutorial.html#classes_list", [
        [ "Types in LatMRG", "tutorial.html#classes_list_types", null ],
        [ "MRG Representation", "tutorial.html#classes_list_mrg", null ],
        [ "Testing and Reducing Lattices", "tutorial.html#classes_list_tests", null ]
      ] ],
      [ "Programming LatMRG", "tutorial.html#program", [
        [ "Working with types", "tutorial.html#program_types", null ],
        [ "Adding new types of generators", "tutorial.html#program_creating", null ],
        [ "Expanding the executable", "tutorial.html#program_expanding", [
          [ "nextGenerator Functions", "tutorial.html#next_gen", null ],
          [ "Accessing the New Function", "tutorial.html#program_exec_mod", null ]
        ] ]
      ] ]
    ] ],
    [ "LatMain", "_lat_main.html", null ],
    [ "SeekMain", "_seek_main.html", [
      [ "SeekMain", "_seek_main.html#seek", [
        [ "Searching method", "_seek_main.html#method_seek", null ],
        [ "The executable programs", "_seek_main.html#executable_seek", null ],
        [ "The data file", "_seek_main.html#data_seek", null ],
        [ "Meaning of the data fields", "_seek_main.html#fields_seek", null ],
        [ "The output files", "_seek_main.html#output_seek", null ]
      ] ]
    ] ],
    [ "Todo List", "todo.html", null ],
    [ "Bibliography", "citelist.html", null ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", null ],
        [ "Functions", "namespacemembers_func.html", null ],
        [ "Variables", "namespacemembers_vars.html", null ],
        [ "Enumerations", "namespacemembers_enum.html", null ],
        [ "Enumerator", "namespacemembers_eval.html", null ]
      ] ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", null ],
        [ "Typedefs", "functions_type.html", null ],
        [ "Enumerations", "functions_enum.html", null ],
        [ "Enumerator", "functions_eval.html", null ],
        [ "Related Functions", "functions_rela.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", null ],
        [ "Functions", "globals_func.html", null ],
        [ "Variables", "globals_vars.html", null ],
        [ "Macros", "globals_defs.html", null ]
      ] ]
    ] ]
  ] ]
];

var NAVTREEINDEX =
[
"_a_w_c_s_w_b_lattice_8h.html",
"class_lat_m_r_g_1_1_m_m_r_g_lattice.html#a31215546f58d221eafd59f5b9461473c",
"classtinyxml2_1_1_mem_pool_t.html#a5fa4fee934a3df2b9e74282244d78390",
"classtinyxml2_1_1_x_m_l_element.html#ae143997e90064ba82326b29a9930ea8f",
"functions_func_g.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';