#include <omp.h>
#include <iostream>
#include <time.h>
#include <fstream>
#include <string>
#include <chrono>#include <ctime>

#define thread_num 4
using namespace std;




struct amino

{
	string name;
	int frequen;
};


void print_codon(string a[],int size )

{ 
	for(int i=0;i<size;i++)
	{
		if(a[i] != "")
			cout<< a[i] <<endl;
	}
}

void print_aminos(amino a[],int size )

{ 
	int i;
	int total=0;



	for( i=0;i<size;i++)
	{
		total += a[i].frequen;
	}



	for(int p=0;p<size;p++)
	{
		cout<<a[p].name<<" frequency is : "<<a[p].frequen<<endl;
		cout<<a[p].name<<" percentage is : "<<(((float)a[p].frequen/(float)total))*100 <<"%"<<endl;
	}

	cout<<"total length of protein : "<<total<<endl;
}


void check_codon_sections(string a[],int size, amino ab[] )
{
	int i;
	int id;
	int bu= size/16;
	omp_set_num_threads(thread_num);
#pragma omp parallel private(i)
	{
#pragma omp sections
		{
#pragma omp section
			{
			
		for( i=0;i<bu;i++)
		{
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}
		}

		#pragma omp section
			{
				
				for( i=bu;i<(2*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
		#pragma omp section
			{
				
				for( i=(2*bu);i<(3*bu);i++)
		{
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	
  #pragma omp section
			{
				
				for( i=(3*bu);i<(4*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(4*bu);i<(5*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(5*bu);i<(6*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(6*bu);i<(7*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(7*bu);i<(8*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(8*bu);i<(9*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(9*bu);i<(10*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(10*bu);i<(11*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(11*bu);i<(12*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(12*bu);i<(13*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(13*bu);i<(14*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(14*bu);i<(15*bu);i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
	 #pragma omp section
			{
				
				for( i=(15*bu);i<size;i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}
 }
	}

}


void check_codon_sequential(string a[],int size, amino ab[] )
{

	for(int i=0;i<size;i++)
	{
		if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
		{
			ab[20].frequen+=1;
			ab[20].name="stopping codon";
			continue;
		}

		if ((a[i]=="ttt") || (a[i]=="ttc"))
		{
			ab[0].frequen+=1;
			ab[0].name="phe";
			continue;
			
		}

		if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
		{
			ab[1].frequen+=1;
			ab[1].name="leu";
			continue;
		}

		if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
		{
			ab[2].frequen+=1;
			ab[2].name="iie";
			continue;
		}

			if ((a[i]=="atg"))
		{
			ab[3].frequen+=1;
			ab[3].name="met";
			continue;
		}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
		{
			ab[4].frequen+=1;
			ab[4].name="val";
			continue;
		}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
		{
			ab[5].frequen+=1;
			ab[5].name="ser";
			continue;
		}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
		{
			ab[6].frequen+=1;
			ab[6].name="pro";
			continue;
		}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
		{
			ab[7].frequen+=1;
			ab[7].name="thr";
			continue;
		}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
		{
			ab[8].frequen+=1;
			ab[8].name="ala";
			continue;
		}

			if ((a[i]=="tat") || (a[i]=="tac"))
		{
			ab[9].frequen+=1;
			ab[9].name="tyr";
			continue;
		}

			if ((a[i]=="cat") || (a[i]=="cac"))
		{
			ab[10].frequen+=1;
			ab[10].name="his";
			continue;
		}

			if ((a[i]=="caa") || (a[i]=="cag"))
		{
			ab[11].frequen+=1;
			ab[11].name="gin";
			continue;
		}
			if ((a[i]=="aat") || (a[i]=="aac"))
		{
			ab[12].frequen+=1;
			ab[12].name="asn";
			continue;
		}

			if ((a[i]=="aaa") || (a[i]=="aag"))
		{
			ab[13].frequen+=1;
			ab[13].name="lys";
			continue;
		}

			if ((a[i]=="gat") || (a[i]=="gac"))
		{
			ab[14].frequen+=1;
			ab[14].name="asp";
			continue;
		}

			if ((a[i]=="gaa") || (a[i]=="gag"))
		{
			ab[15].frequen+=1;
			ab[15].name="glu";
			continue;
		}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
		{
			ab[16].frequen+=1;
			ab[16].name="cys";
			continue;
		}

			if ((a[i]=="tgg"))
		{
			ab[17].frequen+=1;
			ab[17].name="trp";
			continue;
		}

        	if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
		{
			ab[18].frequen+=1;
			ab[18].name="arg";
			continue;
		}

        	if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
		{
			ab[19].frequen+=1;
			ab[19].name="gly";
			continue;
		}

	}


}


void check_codon_parallel(string a[],int size, amino ab[] )
{
	int i;
	int id;
	omp_set_num_threads(thread_num);
#pragma omp parallel
	{
#pragma omp for
		for( i=0;i<size;i++)
		{
			
			if ((a[i]=="taa") || (a[i]=="tag") || (a[i]=="tga"))
			{
#pragma omp critical
				{
					ab[20].frequen+=1;
					ab[20].name="stopping codon";
				}
				continue;
			}

			if ((a[i]=="ttt") || (a[i]=="ttc"))
			{
#pragma omp critical
				{
					ab[0].frequen+=1;
					ab[0].name="phe";
				}
				continue;


			}

			if ((a[i]=="tta") || (a[i]=="ttg") || (a[i]=="ctt") || (a[i]=="ctc") || (a[i]=="cta") || (a[i]=="ctg"))
			{
#pragma omp critical
				{
					ab[1].frequen+=1;
					ab[1].name="leu";
				}
				continue;
			}

			if ((a[i]=="att") || (a[i]=="atc") || (a[i]=="ata"))
			{
#pragma omp critical
				{
					ab[2].frequen+=1;
					ab[2].name="iie";
				}
				continue;
			}

			if ((a[i]=="atg"))
			{
#pragma omp critical
				{
					ab[3].frequen+=1;
					ab[3].name="met";
				}
				continue;
			}


			if ((a[i]=="gtt") || (a[i]=="gtc") || (a[i]=="gta") || (a[i]=="gtg"))
			{
#pragma omp critical
				{
					ab[4].frequen+=1;
					ab[4].name="val";
				}
				continue;
			}

			if ((a[i]=="tct") || (a[i]=="tcc") || (a[i]=="tca") || (a[i]=="tcg") || (a[i]=="agt") || (a[i]=="agc"))
			{
#pragma omp critical
				{
					ab[5].frequen+=1;
					ab[5].name="ser";
				}
				continue;
			}

			if ((a[i]=="cct") || (a[i]=="ccc") || (a[i]=="cca") || (a[i]=="ccg"))
			{
#pragma omp critical
				{
					ab[6].frequen+=1;
					ab[6].name="pro";
				}
				continue;
			}

			if ((a[i]=="act") || (a[i]=="acc") || (a[i]=="aca") || (a[i]=="acg"))
			{
#pragma omp critical
				{
					ab[7].frequen+=1;
					ab[7].name="thr";
				}
				continue;
			}
			if ((a[i]=="gct") || (a[i]=="gcc") || (a[i]=="gca") || (a[i]=="gcg"))
			{
#pragma omp critical
				{
					ab[8].frequen+=1;
					ab[8].name="ala";
				}
				continue;
			}

			if ((a[i]=="tat") || (a[i]=="tac"))
			{
#pragma omp critical
				{
					ab[9].frequen+=1;
					ab[9].name="tyr";
				}
				continue;
			}

			if ((a[i]=="cat") || (a[i]=="cac"))
			{
#pragma omp critical
				{
					ab[10].frequen+=1;
					ab[10].name="his";
				}
				continue;
			}

			if ((a[i]=="caa") || (a[i]=="cag"))
			{
#pragma omp critical
				{
					ab[11].frequen+=1;
					ab[11].name="gin";
				}
				continue;
			}
			if ((a[i]=="aat") || (a[i]=="aac"))
			{
#pragma omp critical
				{
					ab[12].frequen+=1;
					ab[12].name="asn";
				}
				continue;
			}

			if ((a[i]=="aaa") || (a[i]=="aag"))
			{
#pragma omp critical
				{
					ab[13].frequen+=1;
					ab[13].name="lys";
				}
				continue;
			}

			if ((a[i]=="gat") || (a[i]=="gac"))
			{
#pragma omp critical
				{
					ab[14].frequen+=1;
					ab[14].name="asp";
				}
				continue;
			}

			if ((a[i]=="gaa") || (a[i]=="gag"))
			{
#pragma omp critical
				{
					ab[15].frequen+=1;
					ab[15].name="glu";
				}
				continue;
			}

			if ((a[i]=="tgt") || (a[i]=="tgc"))
			{
#pragma omp critical
				{
					ab[16].frequen+=1;
					ab[16].name="cys";
				}
				continue;
			}

			if ((a[i]=="tgg"))
			{
#pragma omp critical
				{
					ab[17].frequen+=1;
					ab[17].name="trp";
				}
				continue;
			}

			if ((a[i]=="cgt") || (a[i]=="cgc") || (a[i]=="cga") || (a[i]=="cgg") || (a[i]=="aga") || (a[i]=="agg"))
			{
#pragma omp critical
				{
					ab[18].frequen+=1;
					ab[18].name="arg";
				}
				continue;
			}

			if ((a[i]=="ggt") || (a[i]=="ggc") || (a[i]=="gga") || (a[i]=="ggg"))
			{
#pragma omp critical
				{
					ab[19].frequen+=1;
					ab[19].name="gly";
				}
				continue;
			}

		}

	}

}

int main()
{ 
	string line;
	string lines [2000];
	string codon;
	string codons [28000];
	string taksir;
	string r;
	amino aminos [21] ;
	int x=0;
	int y=0;
	int k=0;
	int m=1;
	int n=0;
	int j=0;


	for (j=0;j<21;j++)
	{
		aminos[j].frequen=0;
	}
	ifstream myfile ("InputSeq.dat.txt");
	if (myfile.is_open())
	{

		while ( getline (myfile,line) )
		{

			lines[n]=line;

			n++;
		}
	}
	myfile.close();


	for(int z=0;z< (sizeof(lines)/sizeof(*lines));z++)
	{


		int ln=0;
		taksir=lines[z];

		while(m!=taksir.size()){
			r = taksir.substr(y,1);
			y++;
			m++;
			if(r=="a" || r=="c" || r=="t" || r=="g" || r=="u")
			{

				if(ln<3 )
				{

					codon+=r;
					if(ln==2)
					{
						codons[k]=codon;
						k++;
						ln=-1;
						codon="";
					}

					ln++;
				}  
			}    


		} 
		m=0;
		y=0;

	}

	double starts = omp_get_wtime();
	check_codon_sections(codons,(sizeof(codons)/sizeof(*codons)),aminos );
	double ends = omp_get_wtime();
	print_aminos( aminos,21 );
	cout<<"Time of sections in " << thread_num << " threads is :  "<<ends-starts<<endl;

	starts = omp_get_wtime();
	check_codon_parallel(codons,(sizeof(codons)/sizeof(*codons)),aminos );
	 ends = omp_get_wtime();
	cout<<"Time of parallel in " << thread_num << " threads is :  "<<ends-starts<<endl;

	 starts = omp_get_wtime();
	check_codon_sequential(codons,(sizeof(codons)/sizeof(*codons)),aminos );
	 ends = omp_get_wtime();
	cout<<"Time of squential is :  "<<ends-starts<<endl;


	return 0 ;

}



