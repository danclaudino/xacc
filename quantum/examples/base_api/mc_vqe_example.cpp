#include "xacc.hpp"
#include "xacc_service.hpp"
#include "Algorithm.hpp"
#include <fstream>

int main(int argc, char **argv) {
  xacc::Initialize(argc, argv);

  //xacc::set_verbose(true);
  //xacc::logToFile(true);
  //xacc::setLoggingLevel(2);

  const char *data = R"foo(0
Ground state energy: -2263.263771754
Excited state energy: -2263.19429617
Center of mass: -2.820215,16.49446,28.162545
Ground state dipole moment: -2.852215,6.435912,1.520112
Excited state dipole moment: -2.5722,2.3118,0.5654
Transition dipole moment: -0.3489,3.3049,-0.147

1
Ground state energy: -2263.259852808
Excited state energy: -2263.18713945
Center of mass: 1.686348,8.576397,28.655283
Ground state dipole moment: -2.913037,-2.960237,1.929244
Excited state dipole moment: -2.4897,-7.6905,1.5084
Transition dipole moment: -0.3265,3.3791,0.1578

2
Ground state energy: -2263.2610744
Excited state energy: -2263.18942087
Center of mass: 9.257529,3.737151,30.552047
Ground state dipole moment: -1.892731,4.916174,1.860978
Excited state dipole moment: -0.333,3.0429,2.0862
Transition dipole moment: -2.114,2.5778,-0.6786

3
Ground state energy: -2263.273902347
Excited state energy: -2263.2015061
Center of mass: 17.599933,-0.843284,29.936793
Ground state dipole moment: 2.271487,-1.956678,0.350063
Excited state dipole moment: 2.1561,-2.1427,0.2595
Transition dipole moment: 2.3552,-2.2565,-0.3049

4
Ground state energy: -2263.243115856
Excited state energy: -2263.17416282
Center of mass: 27.132018,-0.188324,31.379032
Ground state dipole moment: -2.533506,1.852841,0.872173
Excited state dipole moment: -2.1305,1.5707,0.4719
Transition dipole moment: 3.1982,-0.2978,0.4616

5
Ground state energy: -2263.29019565
Excited state energy: -2263.21759301
Center of mass: 35.520985,2.757106,30.313462
Ground state dipole moment: 8.092243,-1.997972,0.318523
Excited state dipole moment: 0.9679,-2.6629,0.5369
Transition dipole moment: -3.3919,-0.2956,0.2744

6
Ground state energy: -2263.25812333
Excited state energy: -2263.18493523
Center of mass: 41.792566,8.912329,30.264815
Ground state dipole moment: -4.810761,-2.231908,0.920961
Excited state dipole moment: -0.738,0.2717,1.185
Transition dipole moment: -2.7797,-1.617,0.0529

7
Ground state energy: -2263.314629431
Excited state energy: -2263.24431121
Center of mass: 47.917243,15.883589,29.311677
Ground state dipole moment: 1.071047,1.981487,1.090936
Excited state dipole moment: 3.7191,3.137,0.1623
Transition dipole moment: 3.0353,1.5807,-0.6597

8
Ground state energy: -2263.283582177
Excited state energy: -2263.21345739
Center of mass: 47.570923,25.731484,30.279233
Ground state dipole moment: -2.944358,-2.174706,2.37271
Excited state dipole moment: -2.6392,-0.7448,2.3442
Transition dipole moment: -0.8713,-3.2585,-0.1819

9
Ground state energy: -2263.234361844
Excited state energy: -2263.16413403
Center of mass: 46.737516,34.714841,29.373284
Ground state dipole moment: 2.666364,4.833764,1.89874
Excited state dipole moment: 2.448,3.704,1.7472
Transition dipole moment: -0.8853,-3.1402,0.101

10
Ground state energy: -2263.275746815
Excited state energy: -2263.20553908
Center of mass: 40.579436,41.881996,29.191327
Ground state dipole moment: -1.05995,-5.672244,1.66945
Excited state dipole moment: -1.9112,-4.0793,1.6052
Transition dipole moment: -1.5036,3.0183,0.1365

11
Ground state energy: -2263.283549026
Excited state energy: -2263.21341881
Center of mass: 33.708879,48.343194,28.681529
Ground state dipole moment: -1.926097,2.113669,1.619896
Excited state dipole moment: -1.6403,2.0869,1.1687
Transition dipole moment: -1.6069,2.8934,0.1624

12
Ground state energy: -2263.28916608
Excited state energy: -2263.21882588
Center of mass: 25.408293,50.362636,27.936011
Ground state dipole moment: 2.346118,-2.674907,2.082401
Excited state dipole moment: 3.1272,-3.2358,1.8953
Transition dipole moment: -2.9263,1.7835,-0.1373

13
Ground state energy: -2263.272865209
Excited state energy: -2263.20307447
Center of mass: 16.726041,52.22302,27.076359
Ground state dipole moment: -1.554293,3.356681,0.356564
Excited state dipole moment: -2.6373,4.23,-0.088
Transition dipole moment: -2.9889,1.9929,-0.9346

14
Ground state energy: -2263.272137022
Excited state energy: -2263.20127737
Center of mass: 7.920206,47.96269,26.567246
Ground state dipole moment: 4.773509,0.076608,1.331294
Excited state dipole moment: 0.7425,-0.9358,1.4495
Transition dipole moment: -3.1707,-0.8902,0.6489

15
Ground state energy: -2263.290745981
Excited state energy: -2263.2218366
Center of mass: -0.728568,44.103584,25.340287
Ground state dipole moment: -4.162998,0.835023,1.731733
Excited state dipole moment: -4.2095,0.8344,1.4875
Transition dipole moment: -3.3315,0.3335,-1.1616

16
Ground state energy: -2263.274797077
Excited state energy: -2263.20232022
Center of mass: -4.765092,33.601868,27.300097
Ground state dipole moment: 1.092776,2.828762,0.844768
Excited state dipole moment: 1.8467,4.4812,0.2186
Transition dipole moment: 1.3592,3.0472,0.0877

17
Ground state energy: -2263.296869196
Excited state energy: -2263.22028791
Center of mass: -6.033541,24.261661,27.979516
Ground state dipole moment: -5.873466,-3.537204,2.209032
Excited state dipole moment: -2.4219,2.9381,1.1147
Transition dipole moment: 1.5656,2.8158,-0.0976
)foo";

  std::ofstream datafile("datafile.txt", std::ofstream::out);
  datafile << data;
  datafile.close();
  std::string path = "./datafile.txt";
  auto optimizer =
      xacc::getOptimizer("nlopt", {std::make_pair("nlopt-maxeval", 50)});


  std::cout << "TNQVM + ITensor\n\n";
  auto itensor = xacc::getAccelerator("tnqvm",
  {std::make_pair("tnqvm-visitor","itensor-mps"),
                                                std::make_pair("timer",true)});
  auto buffer1 = xacc::qalloc(4);
  auto mc_vqe1 = xacc::getAlgorithm("mc-vqe");
  mc_vqe1->initialize({std::make_pair("accelerator",itensor),
                      std::make_pair("optimizer",optimizer),
                      std::make_pair("data-path",path),
                      std::make_pair("cyclic",true),
                      std::make_pair("log-level",2),
                      std::make_pair("tnqvm-log",false),
                      std::make_pair("nChromophores", 4)});
  mc_vqe1->execute(buffer1);
  auto params = (*buffer1)["opt-params"].as<std::vector<double>>();

  auto buffer2 = xacc::qalloc(4);
  auto p = mc_vqe1->execute(buffer2, params);
 
 for(auto i : p) {
   std::cout << i << "\n";
 }

  /*
  std::cout << "\n\nTNQVM + exaTN\n\n";
  auto exatn =
      xacc::getAccelerator("tnqvm", {std::make_pair("tnqvm-visitor", "exatn"),
                                     std::make_pair("timer", true)});
  auto buffer2 = xacc::qalloc(4);
  auto mc_vqe2 = xacc::getAlgorithm("mc-vqe");
  mc_vqe2->initialize(
      {std::make_pair("accelerator", exatn),
       std::make_pair("optimizer", optimizer),
       std::make_pair("data-path", path), std::make_pair("cyclic", true),
       std::make_pair("log-level", 2), std::make_pair("tnqvm-log", true),
       std::make_pair("nChromophores", 4)});
  mc_vqe2->execute(buffer2);
*/
  ///
  xacc::Finalize();

  return 0;
}
