
=============================================================
Test of correctness and speed of jump-ahead or the MWC. 
The jump ahead is made with the LCG representation, using NTL::ZZ.

=======================================================================
Generator name: baby11, k=1, b=4, m=11, k = 1, a_0 = -1
  m = 11

  Successive states when starting from LCG state y=1, 
  first when running backward, then when running forward. 
    state-y   xx       y   y_{n-k}    xx   next-y 
    mult = 4
       1     [0 1]     1     4     [0 1]     4
       4     [1 1]     4     5     [1 1]     5
       5     [1 2]     5     9     [1 2]     9
       9     [3 0]     9     3     [3 0]     3
       3     [1 0]     3     1     [1 0]     1
       1     [0 1]     1     4     [0 1]     4
       4     [1 1]     4     5     [1 1]     5
       5     [1 2]     5     9     [1 2]     9

    mult = 3
       1     [0 1]     1     4     [0 1]     3
       3     [1 0]     3     1     [1 0]     9
       9     [3 0]     9     3     [3 0]     5
       5     [1 2]     5     9     [1 2]     4
       4     [1 1]     4     5     [1 1]     1
       1     [0 1]     1     4     [0 1]     3
       3     [1 0]     3     1     [1 0]     9
       9     [3 0]     9     3     [3 0]     5

  Perform successive jumps ahead and print MWC state after each jump.
  jumpSize = 3,  numb jumps = 5
  With jumpMWC:
   MWC state =    [3 0]         
   MWC state =    [0 1]         
   MWC state =    [1 2]         
   MWC state =    [1 0]         
   MWC state =    [1 1]         
  With jumpMWCDirect:
   MWC state =    [3 0]         
   MWC state =    [0 1]         
   MWC state =    [1 2]         
   MWC state =    [1 0]         
   MWC state =    [1 1]         

  Perform one large jump ahead of size jumpSize2 = n0 * jumpSize 
  jumpSize2 = 15,  numb jumps = 1, MWC state after the jump: 
   MWC state =    [1 1]         

  Jumps for the MWC state by transforming back and forth to/from the LCG y_{n-k}. 
   Time to make n = 1000000 jumps ahead, jump size = 3,     CPU time:  0.160172
   final LCG state y = 5
   MWC state =    [1 1]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateLagk. 
   Time to make n = 1000000 jumps ahead, jump size = 3,     CPU time:  0.102807
   final LCG state y = 5
   MWC state =    [1 1]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateDirect. 
   Time to make n = 1000000 jumps ahead, jump size = 3,     CPU time:  0.122206
   final LCG state y = 4
   MWC state =    [1 1]         

=======================================================================
Generator name: baby47, k=2, b=4, m=47, k = 2, a_0 = -1
  m = 47

  Successive states when starting from LCG state y=1, 
  first when running backward, then when running forward. 
    state-y   xx       y   y_{n-k}    xx   next-y 
    mult = 4
       1     [0 0 1]     1     16     [0 0 1]     4
       4     [1 0 1]     4     17     [1 0 1]     16
       16     [1 1 1]     16     21     [1 1 1]     17
       17     [1 1 2]     17     37     [1 1 2]     21
       21     [3 1 0]     21     7     [3 1 0]     37
       37     [0 3 1]     37     28     [0 3 1]     7
       7     [2 0 1]     7     18     [2 0 1]     28
       28     [1 2 1]     28     25     [1 2 1]     18
       18     [2 1 0]     18     6     [2 1 0]     25
       25     [0 2 1]     25     24     [0 2 1]     6
       6     [2 0 0]     6     2     [2 0 0]     24
       24     [0 2 0]     24     8     [0 2 0]     2
       2     [0 0 2]     2     32     [0 0 2]     8
       8     [2 0 2]     8     34     [2 0 2]     32
       32     [2 2 2]     32     42     [2 2 2]     34
       34     [3 2 1]     34     27     [3 2 1]     42

    mult = 12
       1     [0 0 1]     1     16     [0 0 1]     12
       12     [0 1 0]     12     4     [0 1 0]     3
       3     [1 0 0]     3     1     [1 0 0]     36
       36     [0 3 0]     36     12     [0 3 0]     9
       9     [3 0 0]     9     3     [3 0 0]     14
       14     [0 1 2]     14     36     [0 1 2]     27
       27     [1 2 0]     27     9     [1 2 0]     42
       42     [2 3 0]     42     14     [2 3 0]     34
       34     [3 2 1]     34     27     [3 2 1]     32
       32     [2 2 2]     32     42     [2 2 2]     8
       8     [2 0 2]     8     34     [2 0 2]     2
       2     [0 0 2]     2     32     [0 0 2]     24
       24     [0 2 0]     24     8     [0 2 0]     6
       6     [2 0 0]     6     2     [2 0 0]     25
       25     [0 2 1]     25     24     [0 2 1]     18
       18     [2 1 0]     18     6     [2 1 0]     28

  Perform successive jumps ahead and print MWC state after each jump.
  jumpSize = 3,  numb jumps = 5
  With jumpMWC:
   MWC state =    [0 1 0]         
   MWC state =    [3 0 0]         
   MWC state =    [2 3 0]         
   MWC state =    [2 0 2]         
   MWC state =    [2 0 0]         
  With jumpMWCDirect:
   MWC state =    [0 1 0]         
   MWC state =    [3 0 0]         
   MWC state =    [2 3 0]         
   MWC state =    [2 0 2]         
   MWC state =    [2 0 0]         

  Perform one large jump ahead of size jumpSize2 = n0 * jumpSize 
  jumpSize2 = 15,  numb jumps = 1, MWC state after the jump: 
   MWC state =    [2 0 0]         

  Jumps for the MWC state by transforming back and forth to/from the LCG y_{n-k}. 
   Time to make n = 1000000 jumps ahead, jump size = 3,     CPU time:  0.288533
   final LCG state y = 25
   MWC state =    [1 2 1]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateLagk. 
   Time to make n = 1000000 jumps ahead, jump size = 3,     CPU time:  0.192714
   final LCG state y = 25
   MWC state =    [1 2 1]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateDirect. 
   Time to make n = 1000000 jumps ahead, jump size = 3,     CPU time:  0.224126
   final LCG state y = 28
   MWC state =    [1 2 1]         

=======================================================================
Generator name: mwc64k2a2, k = 2, a_0 = -1
  m = 2807376099964726257419195024695302342494521706214850559

  Perform successive jumps ahead and print MWC state after each jump.
  jumpSize = 5000,  numb jumps = 5
  With jumpMWC:
   MWC state =    [6341584866925502853 3716243365499047104 198761096088274323]         
   MWC state =    [13123925719800445461 5883763245918026436 400669719351538482]         
   MWC state =    [17600513527753201518 18057907530679789559 533854547054130149]         
   MWC state =    [6656943845141684912 4118482802070459285 205027298425915194]         
   MWC state =    [206203191989042653 15865576399870023030 14230586686974745]         
  With jumpMWCDirect:
   MWC state =    [6341584866925502853 3716243365499047104 198761096088274323]         
   MWC state =    [13123925719800445461 5883763245918026436 400669719351538482]         
   MWC state =    [17600513527753201518 18057907530679789559 533854547054130149]         
   MWC state =    [6656943845141684912 4118482802070459285 205027298425915194]         
   MWC state =    [206203191989042653 15865576399870023030 14230586686974745]         

  Perform one large jump ahead of size jumpSize2 = n0 * jumpSize 
  jumpSize2 = 25000,  numb jumps = 1, MWC state after the jump: 
   MWC state =    [206203191989042653 15865576399870023030 14230586686974745]         

  Jumps for the MWC state by transforming back and forth to/from the LCG y_{n-k}. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.405824
   final LCG state y = 336161551082178912070875276173051554997227137765247797
   MWC state =    [9867919654258248501 9015013476404146328 298601748485540890]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateLagk. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.268823
   final LCG state y = 336161551082178912070875276173051554997227137765247797
   MWC state =    [9867919654258248501 9015013476404146328 298601748485540890]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateDirect. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.383761
   final LCG state y = 1371978343353657022802616054008433934849290887389314973
   MWC state =    [9867919654258248501 9015013476404146328 298601748485540890]         

=======================================================================
Generator name: mwc64k2a2gk, k = 2, a_0 = -886112297149212177
  m = 1762333070936057100131434091229013375173908805022964207

  Perform successive jumps ahead and print MWC state after each jump.
  jumpSize = 5000,  numb jumps = 5
  With jumpMWC:
   MWC state =    [8932889283709195450 15966794459968604284 -728982883275921699]         
   MWC state =    [8417521861719483061 3710729625523886444 -144583557264493350]         
   MWC state =    [7980780708452824421 4455271831306989827 -181373186861008998]         
   MWC state =    [4801368416475705417 16150328739842165065 -756192063582294056]         
   MWC state =    [5696159717963544253 13905904199623452810 -641020630917749301]         
  With jumpMWCDirect:
   MWC state =    [8932889283709195450 15966794459968604284 -728982883275921699]         
   MWC state =    [8417521861719483061 3710729625523886444 -144583557264493350]         
   MWC state =    [7980780708452824421 4455271831306989827 -181373186861008998]         
   MWC state =    [4801368416475705417 16150328739842165065 -756192063582294056]         
   MWC state =    [5696159717963544253 13905904199623452810 -641020630917749301]         

  Perform one large jump ahead of size jumpSize2 = n0 * jumpSize 
  jumpSize2 = 25000,  numb jumps = 1, MWC state after the jump: 
   MWC state =    [5696159717963544253 13905904199623452810 -641020630917749301]         

  Jumps for the MWC state by transforming back and forth to/from the LCG y_{n-k}. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.428772
   final LCG state y = 456568964192057318134526696261433031924768303988148790
   MWC state =    [4744324230349830358 5217770625026484644 -230408633426172203]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateLagk. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.290706
   final LCG state y = 456568964192057318134526696261433031924768303988148790
   MWC state =    [4744324230349830358 5217770625026484644 -230408633426172203]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateDirect. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.376109
   final LCG state y = 498486328660476395322945810197336117601142884219932929
   MWC state =    [4744324230349830358 5217770625026484644 -230408633426172203]         

=======================================================================
Generator name: mwc64k3a2, k = 3, a_0 = -1
  m = 130914317645911944282251326710913006973731472138243469720537068065907539967

  Perform successive jumps ahead and print MWC state after each jump.
  jumpSize = 5000,  numb jumps = 5
  With jumpMWC:
   MWC state =    [3166609257502247023 2757007540422528855 887201392945524416 59344779016031902]         
   MWC state =    [18123259120927090166 16681544166192493692 12197811212902463323 227418483481234465]         
   MWC state =    [6823175858429529396 4070562147729292135 16853403587541239851 85252811337818426]         
   MWC state =    [14250143555036871359 3051132359532863542 5220663423166780762 189428192353414232]         
   MWC state =    [18030454410666045993 12497396368742840972 1940364593687239238 226901002627373717]         
  With jumpMWCDirect:
   MWC state =    [3166609257502247023 2757007540422528855 887201392945524416 59344779016031902]         
   MWC state =    [18123259120927090166 16681544166192493692 12197811212902463323 227418483481234465]         
   MWC state =    [6823175858429529396 4070562147729292135 16853403587541239851 85252811337818426]         
   MWC state =    [14250143555036871359 3051132359532863542 5220663423166780762 189428192353414232]         
   MWC state =    [18030454410666045993 12497396368742840972 1940364593687239238 226901002627373717]         

  Perform one large jump ahead of size jumpSize2 = n0 * jumpSize 
  jumpSize2 = 25000,  numb jumps = 1, MWC state after the jump: 
   MWC state =    [18030454410666045993 12497396368742840972 1940364593687239238 226901002627373717]         

  Jumps for the MWC state by transforming back and forth to/from the LCG y_{n-k}. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.602350
   final LCG state y = 118882376479221297085181944838560870389168662734743299874869431029758037341
   MWC state =    [10976722890207399261 12568111304792052037 11440797799969080269 153097071478145174]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateLagk. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.395021
   final LCG state y = 118882376479221297085181944838560870389168662734743299874869431029758037341
   MWC state =    [10976722890207399261 12568111304792052037 11440797799969080269 153097071478145174]         

  Jumps from LCG state y, with transform to MWC state via LCGtoMWCStateDirect. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.583187
   final LCG state y = 81193961998010715763265562432379517007234367706242830179480740466319304205
   MWC state =    [10976722890207399261 12568111304792052037 11440797799969080269 153097071478145174]         

=======================================================================

Generator name: mcw64k3a3, hardcoded jumps from LCG state y only, no transform to xx. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.051198
   final LCG state y = 118882376479221297085181944838560870389168662734743299874869431029758037341
   MWC state =    [10976722890207399261 12568111304792052037 11440797799969080269 153097071478145174]         

Generator name: mcw64k3a2xx, hardcoded jumps from LCG state y + toMWC, using xx. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.272023
   final LCG state y = 118882376479221297085181944838560870389168662734743299874869431029758037341
   MWC state =    [10976722890207399261 12568111304792052037 11440797799969080269 153097071478145174]         

Generator name: mcw64k3a2, hardcoded jumps from LCG state y, force a_1 = 0. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.250926
   final LCG state y = 118882376479221297085181944838560870389168662734743299874869431029758037341
   MWC state = 10976722890207399261 12568111304792052037 11440797799969080269 153097071478145174

Generator name: mcw64k3a2s, hardcoded jumps from LCG state y + toMWC. 
   Time to make n = 1000000 jumps ahead, jump size = 5000,     CPU time:  0.137499
   final LCG state y = 118882376479221297085181944838560870389168662734743299874869431029758037341
   MWC state = 10976722890207399261 12568111304792052037 11440797799969080269 153097071478145174
