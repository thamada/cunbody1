//Time-stamp: <2008-05-01 09:01:19 hamada>

/*
 * Copyright (C) 2007 
 *      Tsuyoshi Hamada <hamada@progrape.jp>
 *      All rights reserved.
 * This code is released under version 2 of the GNU GPL.
 */



#if (NJ_SHMEM != NPIPE)
DO NOT SUCCESS COMPILATION
#endif


__device__ void
inter_ij(float4 xj,
	 float4 xi,
	 float* mr3i_,
	 float dx_[3])
{
  float dx,dy,dz;
  float r2,r1i,r2i,r3i;
  float mr3i;

  float xj0 = xj.x;
  float xj1 = xj.y;
  float xj2 = xj.z;
  float mj  = xj.w;

  float xi0   = xi.x;
  float xi1   = xi.y;
  float xi2   = xi.z;
  float ieps2 = xi.w;

  dx = xj0 - xi0;
  dy = xj1 - xi1;
  dz = xj2 - xi2;
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

// Extended Chamomile-Scheme
__global__ void
cunbody_kernel(float4* g_xj,
	       float* g_xi,
	       float* g_fi,
	       int ni,
	       int nj)
{
  int i = BID*(blockDim.x)+TID;
  float4 xi;
  float4 OUTPUT = make_float4(0.0, 0.0, 0.0, 0.0);

  // ** Technic ** 
  xi.x = g_xi[i];
  xi.y = g_xi[i+ni];
  xi.z = g_xi[i+ni*2];
  xi.w = g_xi[i+ni*3];

  __shared__ float4 s_xj[NJ_SHMEM];

  for(int j = 0; j<nj; j += 2*NJ_SHMEM){ //------------------ j Loop
    __syncthreads();
    s_xj[TID]          = g_xj[j+TID];          // ** Technic ** 
    s_xj[TID+NJ_SHMEM] = g_xj[j+TID+NJ_SHMEM]; // ** Technic ** 
    __syncthreads();
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
      float mr120, mr121, mr122, mr123, mr124, mr125, mr126, mr127, mr128, mr129;
      float mr130, mr131, mr132, mr133, mr134, mr135, mr136, mr137, mr138, mr139;
      float mr140, mr141, mr142, mr143, mr144, mr145, mr146, mr147, mr148, mr149;
      float mr150, mr151, mr152, mr153, mr154, mr155, mr156, mr157, mr158, mr159;
      float mr160, mr161, mr162, mr163, mr164, mr165, mr166, mr167, mr168, mr169;
      float mr170, mr171, mr172, mr173, mr174, mr175, mr176, mr177, mr178, mr179;
      float mr180, mr181, mr182, mr183, mr184, mr185, mr186, mr187, mr188, mr189;
      float mr190, mr191, mr192, mr193, mr194, mr195, mr196, mr197, mr198, mr199;
      float mr200, mr201, mr202, mr203, mr204, mr205, mr206, mr207, mr208, mr209;
      float mr210, mr211, mr212, mr213, mr214, mr215, mr216, mr217, mr218, mr219;
      float mr220, mr221, mr222, mr223, mr224, mr225, mr226, mr227, mr228, mr229;
      float mr230, mr231, mr232, mr233, mr234, mr235, mr236, mr237, mr238, mr239;
      float mr240, mr241, mr242, mr243, mr244, mr245, mr246, mr247, mr248, mr249;
      float mr250, mr251, mr252, mr253, mr254, mr255;

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
      float dx120[3], dx121[3], dx122[3], dx123[3], dx124[3], dx125[3], dx126[3], dx127[3], dx128[3], dx129[3];
      float dx130[3], dx131[3], dx132[3], dx133[3], dx134[3], dx135[3], dx136[3], dx137[3], dx138[3], dx139[3];
      float dx140[3], dx141[3], dx142[3], dx143[3], dx144[3], dx145[3], dx146[3], dx147[3], dx148[3], dx149[3];
      float dx150[3], dx151[3], dx152[3], dx153[3], dx154[3], dx155[3], dx156[3], dx157[3], dx158[3], dx159[3];
      float dx160[3], dx161[3], dx162[3], dx163[3], dx164[3], dx165[3], dx166[3], dx167[3], dx168[3], dx169[3];
      float dx170[3], dx171[3], dx172[3], dx173[3], dx174[3], dx175[3], dx176[3], dx177[3], dx178[3], dx179[3];
      float dx180[3], dx181[3], dx182[3], dx183[3], dx184[3], dx185[3], dx186[3], dx187[3], dx188[3], dx189[3];
      float dx190[3], dx191[3], dx192[3], dx193[3], dx194[3], dx195[3], dx196[3], dx197[3], dx198[3], dx199[3];

      float dx200[3], dx201[3], dx202[3], dx203[3], dx204[3], dx205[3], dx206[3], dx207[3], dx208[3], dx209[3];
      float dx210[3], dx211[3], dx212[3], dx213[3], dx214[3], dx215[3], dx216[3], dx217[3], dx218[3], dx219[3];
      float dx220[3], dx221[3], dx222[3], dx223[3], dx224[3], dx225[3], dx226[3], dx227[3], dx228[3], dx229[3];
      float dx230[3], dx231[3], dx232[3], dx233[3], dx234[3], dx235[3], dx236[3], dx237[3], dx238[3], dx239[3];
      float dx240[3], dx241[3], dx242[3], dx243[3], dx244[3], dx245[3], dx246[3], dx247[3], dx248[3], dx249[3];
      float dx250[3], dx251[3], dx252[3], dx253[3], dx254[3], dx255[3];

      // ** Technic ** 
      inter_ij(s_xj[0],xi,&mr0,dx0);     inter_ij(s_xj[1],xi,&mr1,dx1);     inter_ij(s_xj[2],xi,&mr2,dx2);     inter_ij(s_xj[3],xi,&mr3,dx3);     inter_ij(s_xj[4],xi,&mr4,dx4);     inter_ij(s_xj[5],xi,&mr5,dx5);     inter_ij(s_xj[6],xi,&mr6,dx6);     inter_ij(s_xj[7],xi,&mr7,dx7);     inter_ij(s_xj[8],xi,&mr8,dx8);     inter_ij(s_xj[9],xi,&mr9,dx9);
      inter_ij(s_xj[10],xi,&mr10,dx10);  inter_ij(s_xj[11],xi,&mr11,dx11);  inter_ij(s_xj[12],xi,&mr12,dx12);  inter_ij(s_xj[13],xi,&mr13,dx13);  inter_ij(s_xj[14],xi,&mr14,dx14);  inter_ij(s_xj[15],xi,&mr15,dx15);  inter_ij(s_xj[16],xi,&mr16,dx16);  inter_ij(s_xj[17],xi,&mr17,dx17);  inter_ij(s_xj[18],xi,&mr18,dx18);  inter_ij(s_xj[19],xi,&mr19,dx19);
      inter_ij(s_xj[20],xi,&mr20,dx20);  inter_ij(s_xj[21],xi,&mr21,dx21);  inter_ij(s_xj[22],xi,&mr22,dx22);  inter_ij(s_xj[23],xi,&mr23,dx23);  inter_ij(s_xj[24],xi,&mr24,dx24);  inter_ij(s_xj[25],xi,&mr25,dx25);  inter_ij(s_xj[26],xi,&mr26,dx26);  inter_ij(s_xj[27],xi,&mr27,dx27);  inter_ij(s_xj[28],xi,&mr28,dx28);  inter_ij(s_xj[29],xi,&mr29,dx29);
      inter_ij(s_xj[30],xi,&mr30,dx30);  inter_ij(s_xj[31],xi,&mr31,dx31);  inter_ij(s_xj[32],xi,&mr32,dx32);  inter_ij(s_xj[33],xi,&mr33,dx33);  inter_ij(s_xj[34],xi,&mr34,dx34);  inter_ij(s_xj[35],xi,&mr35,dx35);  inter_ij(s_xj[36],xi,&mr36,dx36);  inter_ij(s_xj[37],xi,&mr37,dx37);  inter_ij(s_xj[38],xi,&mr38,dx38);  inter_ij(s_xj[39],xi,&mr39,dx39);
      inter_ij(s_xj[40],xi,&mr40,dx40);  inter_ij(s_xj[41],xi,&mr41,dx41);  inter_ij(s_xj[42],xi,&mr42,dx42);  inter_ij(s_xj[43],xi,&mr43,dx43);  inter_ij(s_xj[44],xi,&mr44,dx44);  inter_ij(s_xj[45],xi,&mr45,dx45);  inter_ij(s_xj[46],xi,&mr46,dx46);  inter_ij(s_xj[47],xi,&mr47,dx47);  inter_ij(s_xj[48],xi,&mr48,dx48);  inter_ij(s_xj[49],xi,&mr49,dx49);
      inter_ij(s_xj[50],xi,&mr50,dx50);  inter_ij(s_xj[51],xi,&mr51,dx51);  inter_ij(s_xj[52],xi,&mr52,dx52);  inter_ij(s_xj[53],xi,&mr53,dx53);  inter_ij(s_xj[54],xi,&mr54,dx54);  inter_ij(s_xj[55],xi,&mr55,dx55);  inter_ij(s_xj[56],xi,&mr56,dx56);  inter_ij(s_xj[57],xi,&mr57,dx57);  inter_ij(s_xj[58],xi,&mr58,dx58);  inter_ij(s_xj[59],xi,&mr59,dx59);
      inter_ij(s_xj[60],xi,&mr60,dx60);  inter_ij(s_xj[61],xi,&mr61,dx61);  inter_ij(s_xj[62],xi,&mr62,dx62);  inter_ij(s_xj[63],xi,&mr63,dx63);  inter_ij(s_xj[64],xi,&mr64,dx64);  inter_ij(s_xj[65],xi,&mr65,dx65);  inter_ij(s_xj[66],xi,&mr66,dx66);  inter_ij(s_xj[67],xi,&mr67,dx67);  inter_ij(s_xj[68],xi,&mr68,dx68);  inter_ij(s_xj[69],xi,&mr69,dx69);
      inter_ij(s_xj[70],xi,&mr70,dx70);  inter_ij(s_xj[71],xi,&mr71,dx71);  inter_ij(s_xj[72],xi,&mr72,dx72);  inter_ij(s_xj[73],xi,&mr73,dx73);  inter_ij(s_xj[74],xi,&mr74,dx74);  inter_ij(s_xj[75],xi,&mr75,dx75);  inter_ij(s_xj[76],xi,&mr76,dx76);  inter_ij(s_xj[77],xi,&mr77,dx77);  inter_ij(s_xj[78],xi,&mr78,dx78);  inter_ij(s_xj[79],xi,&mr79,dx79);
      inter_ij(s_xj[80],xi,&mr80,dx80);  inter_ij(s_xj[81],xi,&mr81,dx81);  inter_ij(s_xj[82],xi,&mr82,dx82);  inter_ij(s_xj[83],xi,&mr83,dx83);  inter_ij(s_xj[84],xi,&mr84,dx84);  inter_ij(s_xj[85],xi,&mr85,dx85);  inter_ij(s_xj[86],xi,&mr86,dx86);  inter_ij(s_xj[87],xi,&mr87,dx87);  inter_ij(s_xj[88],xi,&mr88,dx88);  inter_ij(s_xj[89],xi,&mr89,dx89);
      inter_ij(s_xj[90],xi,&mr90,dx90);  inter_ij(s_xj[91],xi,&mr91,dx91);  inter_ij(s_xj[92],xi,&mr92,dx92);  inter_ij(s_xj[93],xi,&mr93,dx93);  inter_ij(s_xj[94],xi,&mr94,dx94);  inter_ij(s_xj[95],xi,&mr95,dx95);  inter_ij(s_xj[96],xi,&mr96,dx96);  inter_ij(s_xj[97],xi,&mr97,dx97);  inter_ij(s_xj[98],xi,&mr98,dx98);  inter_ij(s_xj[99],xi,&mr99,dx99);
      inter_ij(s_xj[100],xi,&mr100,dx100); inter_ij(s_xj[101],xi,&mr101,dx101); inter_ij(s_xj[102],xi,&mr102,dx102); inter_ij(s_xj[103],xi,&mr103,dx103); inter_ij(s_xj[104],xi,&mr104,dx104);
      inter_ij(s_xj[105],xi,&mr105,dx105); inter_ij(s_xj[106],xi,&mr106,dx106); inter_ij(s_xj[107],xi,&mr107,dx107); inter_ij(s_xj[108],xi,&mr108,dx108); inter_ij(s_xj[109],xi,&mr109,dx109);
      inter_ij(s_xj[110],xi,&mr110,dx110); inter_ij(s_xj[111],xi,&mr111,dx111); inter_ij(s_xj[112],xi,&mr112,dx112); inter_ij(s_xj[113],xi,&mr113,dx113); inter_ij(s_xj[114],xi,&mr114,dx114);
      inter_ij(s_xj[115],xi,&mr115,dx115); inter_ij(s_xj[116],xi,&mr116,dx116); inter_ij(s_xj[117],xi,&mr117,dx117); inter_ij(s_xj[118],xi,&mr118,dx118); inter_ij(s_xj[119],xi,&mr119,dx119);
      inter_ij(s_xj[120],xi,&mr120,dx120); inter_ij(s_xj[121],xi,&mr121,dx121); inter_ij(s_xj[122],xi,&mr122,dx122); inter_ij(s_xj[123],xi,&mr123,dx123); inter_ij(s_xj[124],xi,&mr124,dx124);
      inter_ij(s_xj[125],xi,&mr125,dx125); inter_ij(s_xj[126],xi,&mr126,dx126); inter_ij(s_xj[127],xi,&mr127,dx127); inter_ij(s_xj[128],xi,&mr128,dx128); inter_ij(s_xj[129],xi,&mr129,dx129);
      inter_ij(s_xj[130],xi,&mr130,dx130); inter_ij(s_xj[131],xi,&mr131,dx131); inter_ij(s_xj[132],xi,&mr132,dx132); inter_ij(s_xj[133],xi,&mr133,dx133); inter_ij(s_xj[134],xi,&mr134,dx134);
      inter_ij(s_xj[135],xi,&mr135,dx135); inter_ij(s_xj[136],xi,&mr136,dx136); inter_ij(s_xj[137],xi,&mr137,dx137); inter_ij(s_xj[138],xi,&mr138,dx138); inter_ij(s_xj[139],xi,&mr139,dx139);
      inter_ij(s_xj[140],xi,&mr140,dx140); inter_ij(s_xj[141],xi,&mr141,dx141); inter_ij(s_xj[142],xi,&mr142,dx142); inter_ij(s_xj[143],xi,&mr143,dx143); inter_ij(s_xj[144],xi,&mr144,dx144);
      inter_ij(s_xj[145],xi,&mr145,dx145); inter_ij(s_xj[146],xi,&mr146,dx146); inter_ij(s_xj[147],xi,&mr147,dx147); inter_ij(s_xj[148],xi,&mr148,dx148); inter_ij(s_xj[149],xi,&mr149,dx149);
      inter_ij(s_xj[150],xi,&mr150,dx150); inter_ij(s_xj[151],xi,&mr151,dx151); inter_ij(s_xj[152],xi,&mr152,dx152); inter_ij(s_xj[153],xi,&mr153,dx153); inter_ij(s_xj[154],xi,&mr154,dx154);
      inter_ij(s_xj[155],xi,&mr155,dx155); inter_ij(s_xj[156],xi,&mr156,dx156); inter_ij(s_xj[157],xi,&mr157,dx157); inter_ij(s_xj[158],xi,&mr158,dx158); inter_ij(s_xj[159],xi,&mr159,dx159);
      inter_ij(s_xj[160],xi,&mr160,dx160); inter_ij(s_xj[161],xi,&mr161,dx161); inter_ij(s_xj[162],xi,&mr162,dx162); inter_ij(s_xj[163],xi,&mr163,dx163); inter_ij(s_xj[164],xi,&mr164,dx164);
      inter_ij(s_xj[165],xi,&mr165,dx165); inter_ij(s_xj[166],xi,&mr166,dx166); inter_ij(s_xj[167],xi,&mr167,dx167); inter_ij(s_xj[168],xi,&mr168,dx168); inter_ij(s_xj[169],xi,&mr169,dx169);
      inter_ij(s_xj[170],xi,&mr170,dx170); inter_ij(s_xj[171],xi,&mr171,dx171); inter_ij(s_xj[172],xi,&mr172,dx172); inter_ij(s_xj[173],xi,&mr173,dx173); inter_ij(s_xj[174],xi,&mr174,dx174);
      inter_ij(s_xj[175],xi,&mr175,dx175); inter_ij(s_xj[176],xi,&mr176,dx176); inter_ij(s_xj[177],xi,&mr177,dx177); inter_ij(s_xj[178],xi,&mr178,dx178); inter_ij(s_xj[179],xi,&mr179,dx179);
      inter_ij(s_xj[180],xi,&mr180,dx180); inter_ij(s_xj[181],xi,&mr181,dx181); inter_ij(s_xj[182],xi,&mr182,dx182); inter_ij(s_xj[183],xi,&mr183,dx183); inter_ij(s_xj[184],xi,&mr184,dx184);
      inter_ij(s_xj[185],xi,&mr185,dx185); inter_ij(s_xj[186],xi,&mr186,dx186); inter_ij(s_xj[187],xi,&mr187,dx187); inter_ij(s_xj[188],xi,&mr188,dx188); inter_ij(s_xj[189],xi,&mr189,dx189);
      inter_ij(s_xj[190],xi,&mr190,dx190); inter_ij(s_xj[191],xi,&mr191,dx191); inter_ij(s_xj[192],xi,&mr192,dx192); inter_ij(s_xj[193],xi,&mr193,dx193); inter_ij(s_xj[194],xi,&mr194,dx194);
      inter_ij(s_xj[195],xi,&mr195,dx195); inter_ij(s_xj[196],xi,&mr196,dx196); inter_ij(s_xj[197],xi,&mr197,dx197); inter_ij(s_xj[198],xi,&mr198,dx198); inter_ij(s_xj[199],xi,&mr199,dx199);
      inter_ij(s_xj[200],xi,&mr200,dx200); inter_ij(s_xj[201],xi,&mr201,dx201); inter_ij(s_xj[202],xi,&mr202,dx202); inter_ij(s_xj[203],xi,&mr203,dx203); inter_ij(s_xj[204],xi,&mr204,dx204);
      inter_ij(s_xj[205],xi,&mr205,dx205); inter_ij(s_xj[206],xi,&mr206,dx206); inter_ij(s_xj[207],xi,&mr207,dx207); inter_ij(s_xj[208],xi,&mr208,dx208); inter_ij(s_xj[209],xi,&mr209,dx209);
      inter_ij(s_xj[210],xi,&mr210,dx210); inter_ij(s_xj[211],xi,&mr211,dx211); inter_ij(s_xj[212],xi,&mr212,dx212); inter_ij(s_xj[213],xi,&mr213,dx213); inter_ij(s_xj[214],xi,&mr214,dx214);
      inter_ij(s_xj[215],xi,&mr215,dx215); inter_ij(s_xj[216],xi,&mr216,dx216); inter_ij(s_xj[217],xi,&mr217,dx217); inter_ij(s_xj[218],xi,&mr218,dx218); inter_ij(s_xj[219],xi,&mr219,dx219);
      inter_ij(s_xj[220],xi,&mr220,dx220); inter_ij(s_xj[221],xi,&mr221,dx221); inter_ij(s_xj[222],xi,&mr222,dx222); inter_ij(s_xj[223],xi,&mr223,dx223); inter_ij(s_xj[224],xi,&mr224,dx224);
      inter_ij(s_xj[225],xi,&mr225,dx225); inter_ij(s_xj[226],xi,&mr226,dx226); inter_ij(s_xj[227],xi,&mr227,dx227); inter_ij(s_xj[228],xi,&mr228,dx228); inter_ij(s_xj[229],xi,&mr229,dx229);
      inter_ij(s_xj[230],xi,&mr230,dx230); inter_ij(s_xj[231],xi,&mr231,dx231); inter_ij(s_xj[232],xi,&mr232,dx232); inter_ij(s_xj[233],xi,&mr233,dx233); inter_ij(s_xj[234],xi,&mr234,dx234);
      inter_ij(s_xj[235],xi,&mr235,dx235); inter_ij(s_xj[236],xi,&mr236,dx236); inter_ij(s_xj[237],xi,&mr237,dx237); inter_ij(s_xj[238],xi,&mr238,dx238); inter_ij(s_xj[239],xi,&mr239,dx239);
      inter_ij(s_xj[240],xi,&mr240,dx240); inter_ij(s_xj[241],xi,&mr241,dx241); inter_ij(s_xj[242],xi,&mr242,dx242); inter_ij(s_xj[243],xi,&mr243,dx243); inter_ij(s_xj[244],xi,&mr244,dx244);
      inter_ij(s_xj[245],xi,&mr245,dx245); inter_ij(s_xj[246],xi,&mr246,dx246); inter_ij(s_xj[247],xi,&mr247,dx247); inter_ij(s_xj[248],xi,&mr248,dx248); inter_ij(s_xj[249],xi,&mr249,dx249);
      inter_ij(s_xj[250],xi,&mr250,dx250); inter_ij(s_xj[251],xi,&mr251,dx251); inter_ij(s_xj[252],xi,&mr252,dx252); inter_ij(s_xj[253],xi,&mr253,dx253); inter_ij(s_xj[254],xi,&mr254,dx254);
      inter_ij(s_xj[255],xi,&mr255,dx255);

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
	(mr120*dx120[0])+(mr121*dx121[0])+(mr122*dx122[0])+(mr123*dx123[0])+(mr124*dx124[0])+(mr125*dx125[0])+(mr126*dx126[0])+(mr127*dx127[0])+(mr128*dx128[0])+(mr129*dx129[0])+ 
	(mr130*dx130[0])+(mr131*dx131[0])+(mr132*dx132[0])+(mr133*dx133[0])+(mr134*dx134[0])+(mr135*dx135[0])+(mr136*dx136[0])+(mr137*dx137[0])+(mr138*dx138[0])+(mr139*dx139[0])+ 
	(mr140*dx140[0])+(mr141*dx141[0])+(mr142*dx142[0])+(mr143*dx143[0])+(mr144*dx144[0])+(mr145*dx145[0])+(mr146*dx146[0])+(mr147*dx147[0])+(mr148*dx148[0])+(mr149*dx149[0])+ 
	(mr150*dx150[0])+(mr151*dx151[0])+(mr152*dx152[0])+(mr153*dx153[0])+(mr154*dx154[0])+(mr155*dx155[0])+(mr156*dx156[0])+(mr157*dx157[0])+(mr158*dx158[0])+(mr159*dx159[0])+ 
	(mr160*dx160[0])+(mr161*dx161[0])+(mr162*dx162[0])+(mr163*dx163[0])+(mr164*dx164[0])+(mr165*dx165[0])+(mr166*dx166[0])+(mr167*dx167[0])+(mr168*dx168[0])+(mr169*dx169[0])+ 
	(mr170*dx170[0])+(mr171*dx171[0])+(mr172*dx172[0])+(mr173*dx173[0])+(mr174*dx174[0])+(mr175*dx175[0])+(mr176*dx176[0])+(mr177*dx177[0])+(mr178*dx178[0])+(mr179*dx179[0])+ 
	(mr180*dx180[0])+(mr181*dx181[0])+(mr182*dx182[0])+(mr183*dx183[0])+(mr184*dx184[0])+(mr185*dx185[0])+(mr186*dx186[0])+(mr187*dx187[0])+(mr188*dx188[0])+(mr189*dx189[0])+ 
	(mr190*dx190[0])+(mr191*dx191[0])+(mr192*dx192[0])+(mr193*dx193[0])+(mr194*dx194[0])+(mr195*dx195[0])+(mr196*dx196[0])+(mr197*dx197[0])+(mr198*dx198[0])+(mr199*dx199[0])+ 
	(mr200*dx200[0])+(mr201*dx201[0])+(mr202*dx202[0])+(mr203*dx203[0])+(mr204*dx204[0])+(mr205*dx205[0])+(mr206*dx206[0])+(mr207*dx207[0])+(mr208*dx208[0])+(mr209*dx209[0])+ 
	(mr210*dx210[0])+(mr211*dx211[0])+(mr212*dx212[0])+(mr213*dx213[0])+(mr214*dx214[0])+(mr215*dx215[0])+(mr216*dx216[0])+(mr217*dx217[0])+(mr218*dx218[0])+(mr219*dx219[0])+ 
	(mr220*dx220[0])+(mr221*dx221[0])+(mr222*dx222[0])+(mr223*dx223[0])+(mr224*dx224[0])+(mr225*dx225[0])+(mr226*dx226[0])+(mr227*dx227[0])+(mr228*dx228[0])+(mr229*dx229[0])+ 
	(mr230*dx230[0])+(mr231*dx231[0])+(mr232*dx232[0])+(mr233*dx233[0])+(mr234*dx234[0])+(mr235*dx235[0])+(mr236*dx236[0])+(mr237*dx237[0])+(mr238*dx238[0])+(mr239*dx239[0])+ 
	(mr240*dx240[0])+(mr241*dx241[0])+(mr242*dx242[0])+(mr243*dx243[0])+(mr244*dx244[0])+(mr245*dx245[0])+(mr246*dx246[0])+(mr247*dx247[0])+(mr248*dx248[0])+(mr249*dx249[0])+ 
	(mr250*dx250[0])+(mr251*dx251[0])+(mr252*dx252[0])+(mr253*dx253[0])+(mr254*dx254[0])+(mr255*dx255[0]);

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
	(mr120*dx120[1])+(mr121*dx121[1])+(mr122*dx122[1])+(mr123*dx123[1])+(mr124*dx124[1])+(mr125*dx125[1])+(mr126*dx126[1])+(mr127*dx127[1])+(mr128*dx128[1])+(mr129*dx129[1])+ 
	(mr130*dx130[1])+(mr131*dx131[1])+(mr132*dx132[1])+(mr133*dx133[1])+(mr134*dx134[1])+(mr135*dx135[1])+(mr136*dx136[1])+(mr137*dx137[1])+(mr138*dx138[1])+(mr139*dx139[1])+ 
	(mr140*dx140[1])+(mr141*dx141[1])+(mr142*dx142[1])+(mr143*dx143[1])+(mr144*dx144[1])+(mr145*dx145[1])+(mr146*dx146[1])+(mr147*dx147[1])+(mr148*dx148[1])+(mr149*dx149[1])+ 
	(mr150*dx150[1])+(mr151*dx151[1])+(mr152*dx152[1])+(mr153*dx153[1])+(mr154*dx154[1])+(mr155*dx155[1])+(mr156*dx156[1])+(mr157*dx157[1])+(mr158*dx158[1])+(mr159*dx159[1])+ 
	(mr160*dx160[1])+(mr161*dx161[1])+(mr162*dx162[1])+(mr163*dx163[1])+(mr164*dx164[1])+(mr165*dx165[1])+(mr166*dx166[1])+(mr167*dx167[1])+(mr168*dx168[1])+(mr169*dx169[1])+ 
	(mr170*dx170[1])+(mr171*dx171[1])+(mr172*dx172[1])+(mr173*dx173[1])+(mr174*dx174[1])+(mr175*dx175[1])+(mr176*dx176[1])+(mr177*dx177[1])+(mr178*dx178[1])+(mr179*dx179[1])+ 
	(mr180*dx180[1])+(mr181*dx181[1])+(mr182*dx182[1])+(mr183*dx183[1])+(mr184*dx184[1])+(mr185*dx185[1])+(mr186*dx186[1])+(mr187*dx187[1])+(mr188*dx188[1])+(mr189*dx189[1])+ 
	(mr190*dx190[1])+(mr191*dx191[1])+(mr192*dx192[1])+(mr193*dx193[1])+(mr194*dx194[1])+(mr195*dx195[1])+(mr196*dx196[1])+(mr197*dx197[1])+(mr198*dx198[1])+(mr199*dx199[1])+ 
	(mr200*dx200[1])+(mr201*dx201[1])+(mr202*dx202[1])+(mr203*dx203[1])+(mr204*dx204[1])+(mr205*dx205[1])+(mr206*dx206[1])+(mr207*dx207[1])+(mr208*dx208[1])+(mr209*dx209[1])+ 
	(mr210*dx210[1])+(mr211*dx211[1])+(mr212*dx212[1])+(mr213*dx213[1])+(mr214*dx214[1])+(mr215*dx215[1])+(mr216*dx216[1])+(mr217*dx217[1])+(mr218*dx218[1])+(mr219*dx219[1])+ 
	(mr220*dx220[1])+(mr221*dx221[1])+(mr222*dx222[1])+(mr223*dx223[1])+(mr224*dx224[1])+(mr225*dx225[1])+(mr226*dx226[1])+(mr227*dx227[1])+(mr228*dx228[1])+(mr229*dx229[1])+ 
	(mr230*dx230[1])+(mr231*dx231[1])+(mr232*dx232[1])+(mr233*dx233[1])+(mr234*dx234[1])+(mr235*dx235[1])+(mr236*dx236[1])+(mr237*dx237[1])+(mr238*dx238[1])+(mr239*dx239[1])+ 
	(mr240*dx240[1])+(mr241*dx241[1])+(mr242*dx242[1])+(mr243*dx243[1])+(mr244*dx244[1])+(mr245*dx245[1])+(mr246*dx246[1])+(mr247*dx247[1])+(mr248*dx248[1])+(mr249*dx249[1])+ 
	(mr250*dx250[1])+(mr251*dx251[1])+(mr252*dx252[1])+(mr253*dx253[1])+(mr254*dx254[1])+(mr255*dx255[1]);

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
	(mr120*dx120[2])+(mr121*dx121[2])+(mr122*dx122[2])+(mr123*dx123[2])+(mr124*dx124[2])+(mr125*dx125[2])+(mr126*dx126[2])+(mr127*dx127[2])+(mr128*dx128[2])+(mr129*dx129[2])+ 
	(mr130*dx130[2])+(mr131*dx131[2])+(mr132*dx132[2])+(mr133*dx133[2])+(mr134*dx134[2])+(mr135*dx135[2])+(mr136*dx136[2])+(mr137*dx137[2])+(mr138*dx138[2])+(mr139*dx139[2])+ 
	(mr140*dx140[2])+(mr141*dx141[2])+(mr142*dx142[2])+(mr143*dx143[2])+(mr144*dx144[2])+(mr145*dx145[2])+(mr146*dx146[2])+(mr147*dx147[2])+(mr148*dx148[2])+(mr149*dx149[2])+ 
	(mr150*dx150[2])+(mr151*dx151[2])+(mr152*dx152[2])+(mr153*dx153[2])+(mr154*dx154[2])+(mr155*dx155[2])+(mr156*dx156[2])+(mr157*dx157[2])+(mr158*dx158[2])+(mr159*dx159[2])+ 
	(mr160*dx160[2])+(mr161*dx161[2])+(mr162*dx162[2])+(mr163*dx163[2])+(mr164*dx164[2])+(mr165*dx165[2])+(mr166*dx166[2])+(mr167*dx167[2])+(mr168*dx168[2])+(mr169*dx169[2])+ 
	(mr170*dx170[2])+(mr171*dx171[2])+(mr172*dx172[2])+(mr173*dx173[2])+(mr174*dx174[2])+(mr175*dx175[2])+(mr176*dx176[2])+(mr177*dx177[2])+(mr178*dx178[2])+(mr179*dx179[2])+ 
	(mr180*dx180[2])+(mr181*dx181[2])+(mr182*dx182[2])+(mr183*dx183[2])+(mr184*dx184[2])+(mr185*dx185[2])+(mr186*dx186[2])+(mr187*dx187[2])+(mr188*dx188[2])+(mr189*dx189[2])+ 
	(mr190*dx190[2])+(mr191*dx191[2])+(mr192*dx192[2])+(mr193*dx193[2])+(mr194*dx194[2])+(mr195*dx195[2])+(mr196*dx196[2])+(mr197*dx197[2])+(mr198*dx198[2])+(mr199*dx199[2])+ 
	(mr200*dx200[2])+(mr201*dx201[2])+(mr202*dx202[2])+(mr203*dx203[2])+(mr204*dx204[2])+(mr205*dx205[2])+(mr206*dx206[2])+(mr207*dx207[2])+(mr208*dx208[2])+(mr209*dx209[2])+ 
	(mr210*dx210[2])+(mr211*dx211[2])+(mr212*dx212[2])+(mr213*dx213[2])+(mr214*dx214[2])+(mr215*dx215[2])+(mr216*dx216[2])+(mr217*dx217[2])+(mr218*dx218[2])+(mr219*dx219[2])+ 
	(mr220*dx220[2])+(mr221*dx221[2])+(mr222*dx222[2])+(mr223*dx223[2])+(mr224*dx224[2])+(mr225*dx225[2])+(mr226*dx226[2])+(mr227*dx227[2])+(mr228*dx228[2])+(mr229*dx229[2])+ 
	(mr230*dx230[2])+(mr231*dx231[2])+(mr232*dx232[2])+(mr233*dx233[2])+(mr234*dx234[2])+(mr235*dx235[2])+(mr236*dx236[2])+(mr237*dx237[2])+(mr238*dx238[2])+(mr239*dx239[2])+ 
	(mr240*dx240[2])+(mr241*dx241[2])+(mr242*dx242[2])+(mr243*dx243[2])+(mr244*dx244[2])+(mr245*dx245[2])+(mr246*dx246[2])+(mr247*dx247[2])+(mr248*dx248[2])+(mr249*dx249[2])+ 
	(mr250*dx250[2])+(mr251*dx251[2])+(mr252*dx252[2])+(mr253*dx253[2])+(mr254*dx254[2])+(mr255*dx255[2]);

    }

  }//---------------------------------------------------------- j Loop

  // ** Technic ** 
  g_fi[i]      = OUTPUT.x;
  g_fi[i+ni]   = OUTPUT.y;
  g_fi[i+ni*2] = OUTPUT.z;


}


