//Time-stamp: <2008-07-07 20:16:00 hamada>

/*
 * Copyright (C) 2007 
 *      Tsuyoshi Hamada <hamada@progrape.jp>
 *      All rights reserved.
 * This code is released under version 2 of the GNU GPL.
 */

#if (NJ_SHMEM != NPIPE)
DO NOT SUCCESS COMPILATION
#endif

namespace cunbody_kernel_type01{

  __device__ void
  inter(float4 xj,
	   float4 xi,
	   float* mr3i_,
	   float dx_[3])
  {
    float dx,dy,dz;
    float r2,r1i,r2i,r3i;
    float mr3i;

    float mj  = xj.w;
    float ieps2 = xi.w;

    dx = xj.x - xi.x;
    dy = xj.y - xi.y;
    dz = xj.z - xi.z;
    r2 = (dx*dx+ieps2)+ (dy * dy) + (dz * dz); // ** Technic ** 
    r1i = 1/sqrt(r2);
    r2i = r1i * r1i;
    r3i = r2i * r1i;
    mr3i = mj * r3i;

    *mr3i_ = mr3i;
    dx_[0] = dx;
    dx_[1] = dy;
    dx_[2] = dz;
  }


#define TID (threadIdx.x)
#define BID (blockIdx.x)
#define NTH (blockDim.x)

  // Extended Chamomile-Scheme
  __global__ void
  kernel(float4* g_xj,
	 float* g_xi,
	 float* g_fi,
	 int ni,
	 int nj)
  {
    int i = BID*NTH+TID;
    float4 x_i;
    float4 OUTPUT = make_float4(0.0, 0.0, 0.0, 0.0);

    // ** Technic ** 
    x_i.x = g_xi[i];
    x_i.y = g_xi[i+ni];
    x_i.z = g_xi[i+ni*2];
    x_i.w = g_xi[i+ni*3];

    for(int j = 0; j<nj; j += NTH){ //------------------ J-BLOCK Loop
      __shared__ float4 s_xj[NJ_SHMEM];
      __syncthreads();
      s_xj[TID] = g_xj[j+TID];    // ** Technic ** 
      __syncthreads();


#ifdef AC
#undef AC
#endif
#define AC(j,a,b)  inter(s_xj[(j)], x_i, (a), (b));  

      {
	float mr0, mr1, mr2, mr3, mr4, mr5, mr6, mr7, mr8, mr9;
	float mr10, mr11, mr12, mr13, mr14, mr15, mr16, mr17, mr18, mr19;
	float mr20, mr21, mr22, mr23, mr24, mr25, mr26, mr27, mr28, mr29;
	float mr30, mr31, mr32, mr33, mr34, mr35, mr36, mr37, mr38, mr39;
	float mr40, mr41, mr42, mr43, mr44, mr45, mr46, mr47, mr48, mr49;
	float mr50, mr51, mr52, mr53, mr54, mr55, mr56, mr57, mr58, mr59;
	float mr60, mr61, mr62, mr63, mr64, mr65, mr66, mr67, mr68, mr69;
	float mr70, mr71, mr72, mr73, mr74, mr75, mr76, mr77, mr78, mr79;
	float mr80, mr81, mr82, mr83, mr84, mr85, mr86, mr87, mr88, mr89;
	float mr90, mr91, mr92, mr93, mr94, mr95, mr96, mr97, mr98, mr99;
	float mr100, mr101, mr102, mr103, mr104, mr105, mr106, mr107, mr108, mr109;
	float mr110, mr111, mr112, mr113, mr114, mr115, mr116, mr117, mr118, mr119;
	float mr120, mr121, mr122, mr123, mr124, mr125, mr126, mr127;

	float dx0[3], dx1[3], dx2[3], dx3[3], dx4[3], dx5[3], dx6[3], dx7[3], dx8[3], dx9[3];
	float dx10[3], dx11[3], dx12[3], dx13[3], dx14[3], dx15[3], dx16[3], dx17[3], dx18[3], dx19[3];
	float dx20[3], dx21[3], dx22[3], dx23[3], dx24[3], dx25[3], dx26[3], dx27[3], dx28[3], dx29[3];
	float dx30[3], dx31[3], dx32[3], dx33[3], dx34[3], dx35[3], dx36[3], dx37[3], dx38[3], dx39[3];
	float dx40[3], dx41[3], dx42[3], dx43[3], dx44[3], dx45[3], dx46[3], dx47[3], dx48[3], dx49[3];
	float dx50[3], dx51[3], dx52[3], dx53[3], dx54[3], dx55[3], dx56[3], dx57[3], dx58[3], dx59[3];
	float dx60[3], dx61[3], dx62[3], dx63[3], dx64[3], dx65[3], dx66[3], dx67[3], dx68[3], dx69[3];
	float dx70[3], dx71[3], dx72[3], dx73[3], dx74[3], dx75[3], dx76[3], dx77[3], dx78[3], dx79[3];
	float dx80[3], dx81[3], dx82[3], dx83[3], dx84[3], dx85[3], dx86[3], dx87[3], dx88[3], dx89[3];
	float dx90[3], dx91[3], dx92[3], dx93[3], dx94[3], dx95[3], dx96[3], dx97[3], dx98[3], dx99[3];
	float dx100[3], dx101[3], dx102[3], dx103[3], dx104[3], dx105[3], dx106[3], dx107[3], dx108[3], dx109[3];
	float dx110[3], dx111[3], dx112[3], dx113[3], dx114[3], dx115[3], dx116[3], dx117[3], dx118[3], dx119[3];
	float dx120[3], dx121[3], dx122[3], dx123[3], dx124[3], dx125[3], dx126[3], dx127[3];

	// ** Technic ** 
	AC(0,&mr0,dx0);     AC(1,&mr1,dx1);     AC(2,&mr2,dx2);     AC(3,&mr3,dx3);     AC(4,&mr4,dx4);     AC(5,&mr5,dx5);     AC(6,&mr6,dx6);     AC(7,&mr7,dx7);     AC(8,&mr8,dx8);     AC(9,&mr9,dx9);
	AC(10,&mr10,dx10);  AC(11,&mr11,dx11);  AC(12,&mr12,dx12);  AC(13,&mr13,dx13);  AC(14,&mr14,dx14);  AC(15,&mr15,dx15);  AC(16,&mr16,dx16);  AC(17,&mr17,dx17);  AC(18,&mr18,dx18);  AC(19,&mr19,dx19);
	AC(20,&mr20,dx20);  AC(21,&mr21,dx21);  AC(22,&mr22,dx22);  AC(23,&mr23,dx23);  AC(24,&mr24,dx24);  AC(25,&mr25,dx25);  AC(26,&mr26,dx26);  AC(27,&mr27,dx27);  AC(28,&mr28,dx28);  AC(29,&mr29,dx29);
	AC(30,&mr30,dx30);  AC(31,&mr31,dx31);  AC(32,&mr32,dx32);  AC(33,&mr33,dx33);  AC(34,&mr34,dx34);  AC(35,&mr35,dx35);  AC(36,&mr36,dx36);  AC(37,&mr37,dx37);  AC(38,&mr38,dx38);  AC(39,&mr39,dx39);
	AC(40,&mr40,dx40);  AC(41,&mr41,dx41);  AC(42,&mr42,dx42);  AC(43,&mr43,dx43);  AC(44,&mr44,dx44);  AC(45,&mr45,dx45);  AC(46,&mr46,dx46);  AC(47,&mr47,dx47);  AC(48,&mr48,dx48);  AC(49,&mr49,dx49);
	AC(50,&mr50,dx50);  AC(51,&mr51,dx51);  AC(52,&mr52,dx52);  AC(53,&mr53,dx53);  AC(54,&mr54,dx54);  AC(55,&mr55,dx55);  AC(56,&mr56,dx56);  AC(57,&mr57,dx57);  AC(58,&mr58,dx58);  AC(59,&mr59,dx59);
	AC(60,&mr60,dx60);  AC(61,&mr61,dx61);  AC(62,&mr62,dx62);  AC(63,&mr63,dx63);  AC(64,&mr64,dx64);  AC(65,&mr65,dx65);  AC(66,&mr66,dx66);  AC(67,&mr67,dx67);  AC(68,&mr68,dx68);  AC(69,&mr69,dx69);
	AC(70,&mr70,dx70);  AC(71,&mr71,dx71);  AC(72,&mr72,dx72);  AC(73,&mr73,dx73);  AC(74,&mr74,dx74);  AC(75,&mr75,dx75);  AC(76,&mr76,dx76);  AC(77,&mr77,dx77);  AC(78,&mr78,dx78);  AC(79,&mr79,dx79);
	AC(80,&mr80,dx80);  AC(81,&mr81,dx81);  AC(82,&mr82,dx82);  AC(83,&mr83,dx83);  AC(84,&mr84,dx84);  AC(85,&mr85,dx85);  AC(86,&mr86,dx86);  AC(87,&mr87,dx87);  AC(88,&mr88,dx88);  AC(89,&mr89,dx89);
	AC(90,&mr90,dx90);  AC(91,&mr91,dx91);  AC(92,&mr92,dx92);  AC(93,&mr93,dx93);  AC(94,&mr94,dx94);  AC(95,&mr95,dx95);  AC(96,&mr96,dx96);  AC(97,&mr97,dx97);  AC(98,&mr98,dx98);  AC(99,&mr99,dx99);
	AC(100,&mr100,dx100); AC(101,&mr101,dx101); AC(102,&mr102,dx102); AC(103,&mr103,dx103); AC(104,&mr104,dx104); AC(105,&mr105,dx105); AC(106,&mr106,dx106); AC(107,&mr107,dx107); AC(108,&mr108,dx108); AC(109,&mr109,dx109);
	AC(110,&mr110,dx110); AC(111,&mr111,dx111); AC(112,&mr112,dx112); AC(113,&mr113,dx113); AC(114,&mr114,dx114); AC(115,&mr115,dx115); AC(116,&mr116,dx116); AC(117,&mr117,dx117); AC(118,&mr118,dx118); AC(119,&mr119,dx119);
	AC(120,&mr120,dx120); AC(121,&mr121,dx121); AC(122,&mr122,dx122); AC(123,&mr123,dx123); AC(124,&mr124,dx124); AC(125,&mr125,dx125); AC(126,&mr126,dx126); AC(127,&mr127,dx127);

	// ** Technic ** 
	OUTPUT.x += 
	  (mr0*dx0[0])+(mr1*dx1[0])+(mr2*dx2[0])+(mr3*dx3[0])+(mr4*dx4[0])+(mr5*dx5[0])+(mr6*dx6[0])+(mr7*dx7[0])+(mr8*dx8[0])+(mr9*dx9[0])+ 
	  (mr10*dx10[0])+(mr11*dx11[0])+(mr12*dx12[0])+(mr13*dx13[0])+(mr14*dx14[0])+(mr15*dx15[0])+(mr16*dx16[0])+(mr17*dx17[0])+(mr18*dx18[0])+(mr19*dx19[0])+ 
	  (mr20*dx20[0])+(mr21*dx21[0])+(mr22*dx22[0])+(mr23*dx23[0])+(mr24*dx24[0])+(mr25*dx25[0])+(mr26*dx26[0])+(mr27*dx27[0])+(mr28*dx28[0])+(mr29*dx29[0])+ 
	  (mr30*dx30[0])+(mr31*dx31[0])+(mr32*dx32[0])+(mr33*dx33[0])+(mr34*dx34[0])+(mr35*dx35[0])+(mr36*dx36[0])+(mr37*dx37[0])+(mr38*dx38[0])+(mr39*dx39[0])+ 
	  (mr40*dx40[0])+(mr41*dx41[0])+(mr42*dx42[0])+(mr43*dx43[0])+(mr44*dx44[0])+(mr45*dx45[0])+(mr46*dx46[0])+(mr47*dx47[0])+(mr48*dx48[0])+(mr49*dx49[0])+ 
	  (mr50*dx50[0])+(mr51*dx51[0])+(mr52*dx52[0])+(mr53*dx53[0])+(mr54*dx54[0])+(mr55*dx55[0])+(mr56*dx56[0])+(mr57*dx57[0])+(mr58*dx58[0])+(mr59*dx59[0])+ 
	  (mr60*dx60[0])+(mr61*dx61[0])+(mr62*dx62[0])+(mr63*dx63[0])+(mr64*dx64[0])+(mr65*dx65[0])+(mr66*dx66[0])+(mr67*dx67[0])+(mr68*dx68[0])+(mr69*dx69[0])+ 
	  (mr70*dx70[0])+(mr71*dx71[0])+(mr72*dx72[0])+(mr73*dx73[0])+(mr74*dx74[0])+(mr75*dx75[0])+(mr76*dx76[0])+(mr77*dx77[0])+(mr78*dx78[0])+(mr79*dx79[0])+ 
	  (mr80*dx80[0])+(mr81*dx81[0])+(mr82*dx82[0])+(mr83*dx83[0])+(mr84*dx84[0])+(mr85*dx85[0])+(mr86*dx86[0])+(mr87*dx87[0])+(mr88*dx88[0])+(mr89*dx89[0])+ 
	  (mr90*dx90[0])+(mr91*dx91[0])+(mr92*dx92[0])+(mr93*dx93[0])+(mr94*dx94[0])+(mr95*dx95[0])+(mr96*dx96[0])+(mr97*dx97[0])+(mr98*dx98[0])+(mr99*dx99[0])+ 
	  (mr100*dx100[0])+(mr101*dx101[0])+(mr102*dx102[0])+(mr103*dx103[0])+(mr104*dx104[0])+(mr105*dx105[0])+(mr106*dx106[0])+(mr107*dx107[0])+(mr108*dx108[0])+(mr109*dx109[0])+ 
	  (mr110*dx110[0])+(mr111*dx111[0])+(mr112*dx112[0])+(mr113*dx113[0])+(mr114*dx114[0])+(mr115*dx115[0])+(mr116*dx116[0])+(mr117*dx117[0])+(mr118*dx118[0])+(mr119*dx119[0])+ 
	  (mr120*dx120[0])+(mr121*dx121[0])+(mr122*dx122[0])+(mr123*dx123[0])+(mr124*dx124[0])+(mr125*dx125[0])+(mr126*dx126[0])+(mr127*dx127[0]);

	OUTPUT.y += 
	  (mr0*dx0[1])+(mr1*dx1[1])+(mr2*dx2[1])+(mr3*dx3[1])+(mr4*dx4[1])+(mr5*dx5[1])+(mr6*dx6[1])+(mr7*dx7[1])+(mr8*dx8[1])+(mr9*dx9[1])+ 
	  (mr10*dx10[1])+(mr11*dx11[1])+(mr12*dx12[1])+(mr13*dx13[1])+(mr14*dx14[1])+(mr15*dx15[1])+(mr16*dx16[1])+(mr17*dx17[1])+(mr18*dx18[1])+(mr19*dx19[1])+ 
	  (mr20*dx20[1])+(mr21*dx21[1])+(mr22*dx22[1])+(mr23*dx23[1])+(mr24*dx24[1])+(mr25*dx25[1])+(mr26*dx26[1])+(mr27*dx27[1])+(mr28*dx28[1])+(mr29*dx29[1])+ 
	  (mr30*dx30[1])+(mr31*dx31[1])+(mr32*dx32[1])+(mr33*dx33[1])+(mr34*dx34[1])+(mr35*dx35[1])+(mr36*dx36[1])+(mr37*dx37[1])+(mr38*dx38[1])+(mr39*dx39[1])+ 
	  (mr40*dx40[1])+(mr41*dx41[1])+(mr42*dx42[1])+(mr43*dx43[1])+(mr44*dx44[1])+(mr45*dx45[1])+(mr46*dx46[1])+(mr47*dx47[1])+(mr48*dx48[1])+(mr49*dx49[1])+ 
	  (mr50*dx50[1])+(mr51*dx51[1])+(mr52*dx52[1])+(mr53*dx53[1])+(mr54*dx54[1])+(mr55*dx55[1])+(mr56*dx56[1])+(mr57*dx57[1])+(mr58*dx58[1])+(mr59*dx59[1])+ 
	  (mr60*dx60[1])+(mr61*dx61[1])+(mr62*dx62[1])+(mr63*dx63[1])+(mr64*dx64[1])+(mr65*dx65[1])+(mr66*dx66[1])+(mr67*dx67[1])+(mr68*dx68[1])+(mr69*dx69[1])+ 
	  (mr70*dx70[1])+(mr71*dx71[1])+(mr72*dx72[1])+(mr73*dx73[1])+(mr74*dx74[1])+(mr75*dx75[1])+(mr76*dx76[1])+(mr77*dx77[1])+(mr78*dx78[1])+(mr79*dx79[1])+ 
	  (mr80*dx80[1])+(mr81*dx81[1])+(mr82*dx82[1])+(mr83*dx83[1])+(mr84*dx84[1])+(mr85*dx85[1])+(mr86*dx86[1])+(mr87*dx87[1])+(mr88*dx88[1])+(mr89*dx89[1])+ 
	  (mr90*dx90[1])+(mr91*dx91[1])+(mr92*dx92[1])+(mr93*dx93[1])+(mr94*dx94[1])+(mr95*dx95[1])+(mr96*dx96[1])+(mr97*dx97[1])+(mr98*dx98[1])+(mr99*dx99[1])+ 
	  (mr100*dx100[1])+(mr101*dx101[1])+(mr102*dx102[1])+(mr103*dx103[1])+(mr104*dx104[1])+(mr105*dx105[1])+(mr106*dx106[1])+(mr107*dx107[1])+(mr108*dx108[1])+(mr109*dx109[1])+ 
	  (mr110*dx110[1])+(mr111*dx111[1])+(mr112*dx112[1])+(mr113*dx113[1])+(mr114*dx114[1])+(mr115*dx115[1])+(mr116*dx116[1])+(mr117*dx117[1])+(mr118*dx118[1])+(mr119*dx119[1])+ 
	  (mr120*dx120[1])+(mr121*dx121[1])+(mr122*dx122[1])+(mr123*dx123[1])+(mr124*dx124[1])+(mr125*dx125[1])+(mr126*dx126[1])+(mr127*dx127[1]);

	OUTPUT.z += 
	  (mr0*dx0[2])+(mr1*dx1[2])+(mr2*dx2[2])+(mr3*dx3[2])+(mr4*dx4[2])+(mr5*dx5[2])+(mr6*dx6[2])+(mr7*dx7[2])+(mr8*dx8[2])+(mr9*dx9[2])+ 
	  (mr10*dx10[2])+(mr11*dx11[2])+(mr12*dx12[2])+(mr13*dx13[2])+(mr14*dx14[2])+(mr15*dx15[2])+(mr16*dx16[2])+(mr17*dx17[2])+(mr18*dx18[2])+(mr19*dx19[2])+ 
	  (mr20*dx20[2])+(mr21*dx21[2])+(mr22*dx22[2])+(mr23*dx23[2])+(mr24*dx24[2])+(mr25*dx25[2])+(mr26*dx26[2])+(mr27*dx27[2])+(mr28*dx28[2])+(mr29*dx29[2])+ 
	  (mr30*dx30[2])+(mr31*dx31[2])+(mr32*dx32[2])+(mr33*dx33[2])+(mr34*dx34[2])+(mr35*dx35[2])+(mr36*dx36[2])+(mr37*dx37[2])+(mr38*dx38[2])+(mr39*dx39[2])+ 
	  (mr40*dx40[2])+(mr41*dx41[2])+(mr42*dx42[2])+(mr43*dx43[2])+(mr44*dx44[2])+(mr45*dx45[2])+(mr46*dx46[2])+(mr47*dx47[2])+(mr48*dx48[2])+(mr49*dx49[2])+ 
	  (mr50*dx50[2])+(mr51*dx51[2])+(mr52*dx52[2])+(mr53*dx53[2])+(mr54*dx54[2])+(mr55*dx55[2])+(mr56*dx56[2])+(mr57*dx57[2])+(mr58*dx58[2])+(mr59*dx59[2])+ 
	  (mr60*dx60[2])+(mr61*dx61[2])+(mr62*dx62[2])+(mr63*dx63[2])+(mr64*dx64[2])+(mr65*dx65[2])+(mr66*dx66[2])+(mr67*dx67[2])+(mr68*dx68[2])+(mr69*dx69[2])+ 
	  (mr70*dx70[2])+(mr71*dx71[2])+(mr72*dx72[2])+(mr73*dx73[2])+(mr74*dx74[2])+(mr75*dx75[2])+(mr76*dx76[2])+(mr77*dx77[2])+(mr78*dx78[2])+(mr79*dx79[2])+ 
	  (mr80*dx80[2])+(mr81*dx81[2])+(mr82*dx82[2])+(mr83*dx83[2])+(mr84*dx84[2])+(mr85*dx85[2])+(mr86*dx86[2])+(mr87*dx87[2])+(mr88*dx88[2])+(mr89*dx89[2])+ 
	  (mr90*dx90[2])+(mr91*dx91[2])+(mr92*dx92[2])+(mr93*dx93[2])+(mr94*dx94[2])+(mr95*dx95[2])+(mr96*dx96[2])+(mr97*dx97[2])+(mr98*dx98[2])+(mr99*dx99[2])+ 
	  (mr100*dx100[2])+(mr101*dx101[2])+(mr102*dx102[2])+(mr103*dx103[2])+(mr104*dx104[2])+(mr105*dx105[2])+(mr106*dx106[2])+(mr107*dx107[2])+(mr108*dx108[2])+(mr109*dx109[2])+ 
	  (mr110*dx110[2])+(mr111*dx111[2])+(mr112*dx112[2])+(mr113*dx113[2])+(mr114*dx114[2])+(mr115*dx115[2])+(mr116*dx116[2])+(mr117*dx117[2])+(mr118*dx118[2])+(mr119*dx119[2])+ 
	  (mr120*dx120[2])+(mr121*dx121[2])+(mr122*dx122[2])+(mr123*dx123[2])+(mr124*dx124[2])+(mr125*dx125[2])+(mr126*dx126[2])+(mr127*dx127[2]);
      }

    }//---------------------------------------------------------- J-BLOCK Loop

    // ** Technic ** 
    g_fi[i]      = OUTPUT.x;
    g_fi[i+ni]   = OUTPUT.y;
    g_fi[i+ni*2] = OUTPUT.z;

  }


};
