/*
* MC.cpp
* Simulation of 3D Rods By GCMC with SUS, No energy
* Author: Yuding Ai
* Date: 2015.12.01
* *************************** MC implementation ********************************
* This simulation follows Successive Umbrella sampling on order parameter S.
* ******************************************************************************
*
*/


#include "MC.h"

MC::MC(long int ST, int LEN, int N0, int N1, int N2, double Z)
{
	VRodlist.clear(); // the list storage the Vertical Rods;
    HRodlist.clear(); // the list storage the Horizantal Rods;
    URodlist.clear(); // the list storage the Up Rods;
	n0 = N0;          // length of the box
	n1 = N1;          // weight of the box
	n2 = N2;          // hight  of the box
	length = LEN;     // length of rod
	step = ST;        // MC steps
	z = Z;            // activity z = exp(beta*miu)
	nh=nv=nu=dh=dv=du=ah=av=au=0;            // initialize all my number count into 0; ie. nh stands for number of horizontal rods
	                                         // dh stands for the times of deletion for for horizontal rods and ah is for addition.
}

const vector<HR>& MC::getVRodlist() const
{
	return VRodlist;
}

const vector<HR>& MC::getHRodlist() const
{
	return HRodlist;
}

const vector<HR>& MC::getURodlist() const
{
	return URodlist;
}


void MC::Add(Cells &s,double &prob,double &proba)
{
	int x,y,z,o; // pick a random position and orientation for the HR to be added;
	x = rand()%n0; // x range[0,n0-1]
	y = rand()%n1; // y range[0,n1-1]
	z = rand()%n2; // z range[0,n2-1]
	o = rand()%3  ; // 0 range {0,1,2}

	if(s.getSquare(x,y,z).isEmpty()) // if it's open, try to do Addition;
	{
		HR rod(x,y,z,length,o);

		if(o == 0)
		{
		//======================== Vertical ===============================
			// the vertical case
			int counter = 0;

			for (int j = 0; j < length-1; j++)
			{
				// check if the vertical space is wide open
				if(s.getSquare(x,(y+j+1)%n1,z).isOccupied())
				{
					counter++;
					break;
				}						
			}
			// cout <<"counter = "<< counter<< endl;
			if (counter == 0)
			{
				if(prob<=proba)
				{
					// Do addition;
					// push the new rod into the VRodlist;
					VRodlist.push_back(rod);
					av++;
					nv++;// accumulate the # of ver rod;
					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{	
						s.getSquare(x,(y+i)%n1,z).setStatus(1);
					}					
				}		
			}									
		}

		else if(o == 1)
		{
        //======================= Horizontal  ============================
			int counter = 0;
			for (int j = 0; j< length-1 ; j++)
			{
				// check if the horizontal space is wide open
				if(s.getSquare((x+1+j)%n0,y,z).isOccupied())
				{
					counter++;
					break;
				}							
			}
			if (counter == 0)
			{
				if(prob<= proba)
				{
					//Do addition;
					//push the new rod into the HRodlist;
					HRodlist.push_back(rod);
					ah++;
					nh++;// accumulate the # of hor rod;

					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{
						s.getSquare((x+i)%n0,y,z).setStatus(1);
					}
				}
			}
		}
		else 
		{
        //======================= Up  ============================
			int counter = 0;
			for (int j = 0; j< length-1 ; j++)
			{
				// check if the horizontal space is wide open
				if(s.getSquare(x,y,(z+j+1)%n2).isOccupied())
				{
					counter++;
					break;
				}							
			}
			if (counter == 0)
			{
				if(prob<= proba)
				{
					//Do addition;
					//push the new rod into the HRodlist;
					URodlist.push_back(rod);
					au++;
					nu++;// accumulate the # of hor rod;

					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{
						s.getSquare(x,y,(z+i)%n2).setStatus(1);
					}
				}
			}							
		}
    }
}

void MC::Del(Cells &s,double &prob,double &probd, double &size)
{

	if(nv + nh + nu >0)// make sure there are Vertical rod;
	{
		int indx; // pick a random index from the Rodlist;
		indx = rand()%int(nv+nh+nu);

		//remove Rodlist[indx];
		int x,y,z;// the position of the target on the cells;

		if(prob <= probd)
		{					
			if(indx < nv)
			{
				x = VRodlist[indx].getX();
				y = VRodlist[indx].getY();
				z = VRodlist[indx].getZ();	
				// --------------------- it's a vertical rod -----------------------

				for(int i = 0; i<VRodlist[indx].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare(x,(y+i)%n1,z).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				VRodlist.erase(VRodlist.begin() + indx);
				nv--;// substract the # of ver rod;
				dv++;
			}
			
			else if (indx < nv + nh)
			{
				x = HRodlist[indx - nv].getX();
				y = HRodlist[indx - nv].getY();
				z = HRodlist[indx - nv].getZ();	
				// --------------------- it's a horizontal rod -----------------------

				for(int i = 0; i<HRodlist[indx-nv].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare((x+i)%n0,y,z).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				HRodlist.erase(HRodlist.begin() + indx - nv);
				nh--;// substract the # of ver rod;
				dh++;
			}
			else
			{
				x = URodlist[indx - nv - nh].getX();
				y = URodlist[indx - nv - nh].getY();
				z = URodlist[indx - nv - nh].getZ();	
				// --------------------- it's a up rod -----------------------
				for(int i = 0; i<URodlist[indx-nv-nh].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare(x,y,(z+i)%n2).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				URodlist.erase(URodlist.begin()+indx-nv-nh);
				nu--;// substract the # of ver rod;
				du++;
			}
		}
	}	
}

/*
***************************** The idea of weight on S **************************
* Dividing the order parameter axis into equally sized intervals of width DeltaS
* S = (rho_3-(rho_1+rho_2)/2)/(rho_1+rho_2+rho_3) ; rho_i = ni/V
*
* HERE: We pick rho_3 = rho_u; So the nematic phase has to be up
*
* NOTICED The overlap of the window! As:
* the 1st nematic order parameter is confined to 0 <= S < DeltaS, 
* then 2nd: DeltaS/2 <= S < 3/2*(DeltaS)
* ...  ith: (i-1)*DeltaS/2 <= S < (i+1)*DeltaS/2 
*
* According to R.L.C.Vink's paper, Delta S = 0.001 to 0.002 gives good result
*********************************************************************************
*/	

array<double,100000> MC::MCSUS()
{
	Cells s(n0,n1,n2,EMPTY,length);  //  setting the lattice;
	//================= declare instance variables ============================= 
	stringstream sh,ph;
	sh.precision(20);
	ph.precision(20);
	double addordel;           // the prob to decide either add or del;
	double proba,probd; // the acceptance prob of addition; proba = min(1.0,aaccp);
	// double probdh, probdv, probdu; // the acceptance prob of deletion; probd = min(1.0,daccp);
	double prob; // the prob to decide either accept add/del;
	double aaccp,daccp; 
	double tho; // the density 
	double size;
	double V = double(n0*n1*n2);    // the total lattice size
	double S = -1;

	srand(time(NULL));
	long int i = 0; // counter for each window
	double w = 1.0; // a counter that keep track of the index of window
	double fu,fl; // occurrence counter
	double DeltaS = 0.001; // Declare the DeltaS Here!
	double N_window = 0.5*1/(0.5*DeltaS);
    array<double,100000> WF; // NOTICED THAT THE # OF WEIGHTS IS NOT THE TOTAL NUMBER ANYMORE!!!!
	array<double,100000> P_w; // 
	array<double,100000> PH_w;
	array<double,100000> PL_w;

    // ============= Do a GCMC to reach the initial state where S is in [0,deltaS] ==========
	while (S < 0 || S>DeltaS)
	{
		S = (nu-(nv+nh)*0.5)/(nu+nv+nh);
		addordel = rand()%2 ; 
		size = av+ah+au-dv-dh-du;

		// *****************define the probabilities ***********************************
		prob = ((double) rand() / (RAND_MAX)); 
		tho = double(length*size)/double(n0*n1*n2);

		aaccp = z*double(n0*n1*n2)/(double(nh+nv+nu+1.0)*double(length));
		daccp = (double(nh+nv+nu)*double(length))/(z*double(n0*n1*n2));

		proba = min(1.0,aaccp);
		probd = min(1.0,daccp);
        // ===========================Addition ===================================
		if(addordel == 0) 
		{
			Add(s,prob,proba);
		}
		// ============================Deletion=============================
		else
		{
			if (size != 0) // make sure there are rods to be del;
			{
				Del(s,prob,probd,size);
			}
		}
		// cout << S << endl;			
	}
	// cout << "LOL" << endl;
	//================================Start my SUS_MC simulation=================================
    // Weights WF = [(s=0.0005),(s=0.001),(s=0.0015),(s=0.002) ... , (s=0.5)] ---> say 0.5/0.0005 = 1000 components
    // Windows = [{(s=0.0005),(s=0.001)}; {(s=0.001),(s=0.0015)};...;{(s=0.4995),(s=0.5)}] ---> 0.4995/0.0005 = 999 windows
    // !! the while loop start from a state that when S is in [0,deltaS] !! NOT EMPTY!!!
	//===========================================================================================

	// ================================== THE SUS on S simulation starts ==================================
	while (w <= N_window)// the while loop terminates when finish simulating the last window: window[V/K]; !!!for rods: window[0.8*V/K]
	{
		i = 0;  // initialize my step counter; 
		fl = fu = 0; // initialize my occurrence counter;

		while(i < step )
		// Simulation for each window with "step" amount of step		
		{
			i++;
			// generate a random probability to decide either add or del;
			addordel = rand()%2;
		    double size = nv+nh+nu;
		    double N = int(size);


			prob = ((double) rand() / (RAND_MAX)); 

	        // ===========================Addition ===================================
			if(addordel == 0) 
			{
				Add(s,prob,proba);
			}
			// ============================Deletion=============================
			else 
			{
				Del(s,prob,probd,size);	
			}	

			S = (nu-(nv+nh)*0.5)/(nu+nv+nh);

			// Now check S, only accept move when S is inside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
			// if ((w-1.0)*DeltaS/2.0<= S && S <= (w+1.0)*DeltaS/2.0)
			// {
			if ((w/2.0)*DeltaS <= S && S <= (w+1.0)*DeltaS/2.0) // update the fu after each step.
			{
				fu++; // if at the upper window, update fu
				PH_w[N]++; // update the distribution/histogram of PH;
				cout << "!!" <<nv<<"   "<<nh<<"   "<<nu<<"     "<<S<<endl;
			}
			else if ((w-1.0)*DeltaS/2.0 <= S && S < (w/2.0)*DeltaS)  // update the fl after each step.
			{
				fl++;//if at the lower window, update fl
				PL_w[N]++; // update the distribution/histogram of PL;
				cout << "??"<<endl;
			}	
			// }
			else
			{
				
				if(addordel == 0) 
				{
					Add_redo();
				}
				// ============================Deletion=============================
				else 
				{
					Del_redo();
				}	
			} // the case that S is outside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
		}
		// =======================  if fu and fl != 0 Update the upper window side ================================
        if (fu!=0 && fl != 0)
        {
		    WF[w] = WF[w] + log(fu/fl);
	        // "linearly extrapolate" for WF[w+1] by using W[w] and WF[w-1]
	        WF[w+1] = WF[w];
	        // WF[w+1] = 2*WF[w] - WF[w-1];

	        cout << nu<<"  "<<nh <<" " <<nv <<"  "<<(nu-(nv+nh)*0.5)/(nu+(nv+nh)*0.5)<<endl;
			// ======================= Print out the data into terminal =============================================		
			cout <<"Window: "<< w <<" : "<<"W("<<w<<" : lower) = "<< WF[w-1]<<" "<<"W("<<w<<" : Upper) = "<< WF[w] << endl;
			// initial config determine the intial value of fu and fl
			// Store the distribution of PH_w(N) and PL_w(N)
		    w++; // switch into the next window
        }
        // else do nothing and reset fu and fl and repeat the window simulation
	}

	//accumulate the distribution PH_w(N) and PL_w(N) into the total distribution of P_w(N);
	// for(int k = 0; k< V/int(length) + 1;k++)
	// {
	// 	P_w[k] = PH_w[k] + PL_w[k+1]; 
	// 	ph<<P_w[k]<<endl;
	// }

	for(int i = 0; i< N_window+1; i++) 
	{
		sh<<WF[i]<<endl;
	}

	ofstream myfile ("SUSWeight_function.txt");
	myfile.precision(20);
	string data2 = sh.str();
	// string data3 = ph.str();
	myfile << data2;
	// <<"   "<<data3;
	myfile.close();

	return WF;  
}





