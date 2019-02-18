from numpy import array , zeros
# These are Equation of state parameters for Calculation of state compression factor/Composition analysis
a = array([.153832600, 1.341953000, -2.998583000, -.048312280, .375796500, -1.589575000, -.053588470, .886594630,
      -.710237040, -1.471722000, 1.321850350, -.786659250, 2.291290*10**(-9), .157672400, -.436386400, -0.044081590,
      -.003433888, .032059050, .024873550, .073322790, -.001600573, .642470600, -.416260100, -.066899570,
      .279179500, -.696605100, -.002860589, -.008098836, 3.150547000, .007224479, -.705752900, .534979200,
      -.079314910, -1.418465000, -5.99905*10**(-17), .105840200, .034317290, -.007022847, .024955870, .042968180,
      .746545300, -.291961300, 7.294616000, -9.936757000, -.005399808, -.243256700, .049870160, .003733797, 1.874951000,
      .002168144, -.658716400, .000205518, .009776195, -.020487080, .015573220, .006862415, -.001226752, .002850908 ])

b = array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
      4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 6, 6, 7, 7, 8, 8, 8, 9, 9 ])

c = array(12 * [0] + 6 * [1] + 2 * [0] + 7 * [1] + [0] + 9 * [1] + 2 * [0] + 5 * [1] + [0] + 4 * [1] + [0, 1, 0] + 6 * [1])

k = array([0, 3, 2, 2, 2, 4, 4, 0, 0, 2, 2, 2, 4, 4, 4, 4, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 0, 0, 2, 2, 2, 4, 4, 0, 2, 2, 4, 4, 0, 2, 0, 2, 1, 2, 2, 2, 2 ])

u = array([0.0, 0.5, 1.0, 3.5, -.5, 4.5, .5, 7.5, 9.5, 6.0, 12.0, 12.5, -6.0, 20, 3.0, 2.0, 2.0, 11.0, -.5, .5, 0.0, 4.0, 6.0, 21.0, 23.0, 22.0, -1.0,
      -.5, 7.0, -1.0, 6.0, 4.0, 1.0, 9.0, -13.0, 21.0, 8.0, -.5, 0.0, 2.0, 7.0, 9.0, 22.0, 23.0, 1.0, 9.0, 3.0, 8.0, 23.0, 1.5, 5.0, -.5, 4.0, 7.0, 3.0, 0.0, 1.0, 0.0 ])

g = array(4 * [0] + [1, 1] + 18 * [0] + [1, 0, 0, 0, 1, 0, 0, 1, 1, 1] + 16 * [0] + [1, 0, 0, 1, 0, 1, 0, 0])

q = zeros(58) ; q[6]= 1 ; q[15]= 1 ; q[25]= 1 ; q[27]= 1 ; q[36]= 1 ; q[41]= 1 ; q[46]= 1 ; q[48]= 1 ; q[51]= 1 ; q[57]= 1 

f = zeros(58) ; f[12] = 1 ; f[26] = 1 ; f[29] = 1 ; f[34] = 1 

s = zeros(58) ; s[7] = 1 ; s[8] = 1

w = zeros(58) ; w[9] = 1 ; w[10] = 1 ; w[11] = 1 