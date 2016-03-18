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
	o = rand()%2  ; // 0 range {0,1,2}

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



void MC::AddSUS(Cells &s,double &prob,double &proba,double &w,double &DeltaS)
{
	int x,y,z,o; // pick a random position and orientation for the HR to be added;
	x = rand()%n0; // x range[0,n0-1]
	y = rand()%n1; // y range[0,n1-1]
	z = rand()%n2; // z range[0,n2-1]
	o = rand()%3  ; // 0 range {0,1,2}
	double S = 0;

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
					S = (nu-(nv+1+nh)*0.5)/(nu+nv+1+nh);

					// Now check S, only accept move when S is inside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
					if (-0.5+(w-1.0)*DeltaS/2.0 <= S && S <= -0.5+(w+1.0)*DeltaS/2.0) // update the fu after each step.
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
					    return ;	
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
					S = (nu-(nv+1+nh)*0.5)/(nu+nv+1+nh);

					// Now check S, only accept move when S is inside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
					if (-0.5+(w-1.0)*DeltaS/2.0 <= S && S <= -0.5+(w+1.0)*DeltaS/2.0) 
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
					    return ;
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
					S = (nu+1-(nv+nh)*0.5)/(nu+nv+1+nh);

					// Now check S, only accept move when S is inside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
					if (-0.5+(w-1.0)*DeltaS/2.0 <= S && S <= -0.5+(w+1.0)*DeltaS/2.0) 
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
					    return ;
					}
				}
			}							
		}
    }
    return ;
}


void MC::DelSUS(Cells &s,double &prob,double &probd, double &size,double&w,double&DeltaS)
{		
	int x,y,z,indx;// the position of the target on the cells;
	x=y=z=indx=0;
	double S = 0;

	if(nv + nh + nu >0)// make sure there are Vertical rod;
	{	
		if(prob <= probd)
		{
			indx = rand()%int(nv+nh+nu);; // pick a random index from the Rodlist;
			if(indx < nv && nv !=0)
			{
				x = VRodlist[indx].getX();
				y = VRodlist[indx].getY();
				z = VRodlist[indx].getZ();	

				// --------------------- it's a vertical rod -----------------------
				S = (nu-(nv-1+nh)*0.5)/(nu+nv-1+nh);

				// Now check S, only accept move when S is inside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
				if (-0.5+(w-1.0)*DeltaS/2.0 <= S && S <= -0.5+(w+1.0)*DeltaS/2.0) 
				{
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
			}
			
			else if (indx < nv + nh && nh !=0)
			{
				x = HRodlist[indx - nv].getX();
				y = HRodlist[indx - nv].getY();
				z = HRodlist[indx - nv].getZ();	

				// --------------------- it's a horizontal rod -----------------------
				S = (nu-(nv-1+nh)*0.5)/(nu+nv-1+nh);

				// Now check S, only accept move when S is inside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
				if (-0.5+(w-1.0)*DeltaS/2.0 <= S && S <= -0.5+(w+1.0)*DeltaS/2.0) 
				{
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
			}
			else if(nu !=0)
			{
				x = URodlist[indx - nv - nh].getX();
				y = URodlist[indx - nv - nh].getY();
				z = URodlist[indx - nv - nh].getZ();	

				// --------------------- it's a up rod -----------------------
				S = (nu-1-(nv+nh)*0.5)/(nu+nv-1+nh);

				// Now check S, only accept move when S is inside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
				if (-0.5+(w-1.0)*DeltaS/2.0 <= S && S <= -0.5+(w+1.0)*DeltaS/2.0) 
				{
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
			return ;					
		}
		return ; // if delesion failed, return -1,-1,-1,-1
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
* At the end, for Window i: we should have P_wi(N) = P_H_wi(N)+P_L_wi+1(N)
* and we obtained the P(N) by P(N) = C*sum_w1->wmax{P_wi(N)*exp{WF(i)}}
*
* According to R.L.C.Vink's paper, Delta S = 0.001 to 0.002 gives good result
*********************************************************************************
*/	

// change to return a tuple:<WF,P_w(a two dimensional array that stores pi(N))>

// array<double,2000> MC::MCSUS()
void MC::MCSUS()
// std::tuple<array<double,2000>,vector<vector<double>>> MC::MCSUS()
{
	Cells s(n0,n1,n2,EMPTY,length);  //  setting the lattice;
	//================= declare instance variables ============================= 
	stringstream sh,ph,data_s,extrap;
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
	long int j = 0;
	double w = 1.0; // a counter that keep track of the index of window
	double fu,fl; // occurrence counter
	double DeltaS = 0.001; // Declare the DeltaS Here!
	double N_window = 0.8*1.3/(0.5*DeltaS);
    array<double,3000> WF; // NOTICED THAT THE # OF WEIGHTS IS NOT THE TOTAL NUMBER ANYMORE!!!!   --> 1D array
    vector<double> P_N;
    // vector<double> P_N_Z2;
    vector<double> CN; // an array of normalization constant for P_w
    // P_ext.resize(n0*n1*n2/length);
    P_N.resize(n0*n1*n2/length);
    // P_N_Z2.resize(n0*n1*n2/length);
    CN.resize(n0*n1*n2/length);
    // -----------------------The rest arrays are 2D arrays ------------------------------//
	vector<vector<double>> P_w; // an array that stores the histogram/distribution of ave # of N
	                            // where P_w^{i}(N) = PH_w^{i}(N) + PL_w^{i+1}(N);
	// vector<vector<double>> P_w_Z2;
	vector<vector<double>> PH_w; // an array that stores the histogram/distribution of Higher window in terms of # of N
	vector<vector<double>> PL_w; // an array that stores the histogram/distribution of Lower window in terms of # of N
	PH_w.resize(3000);
	PL_w.resize(3000);
	P_w.resize(3000);
	// P_w_Z2.resize(3000);

	for(int i =0; i<3000; i++)
	{
		PH_w[i].resize(n0*n1*n2/length);
		PL_w[i].resize(n0*n1*n2/length);
		P_w[i].resize(n0*n1*n2/length);
		// P_w_Z2[i].resize(n0*n1*n2/length);
		// cout<< "LOL: "<<i<<endl;
	}

    // ============= Do a GCMC to reach a initial state where S is -0.5==========
	
	while(j<step)
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
		j++;

		// cout << S << endl;			
	}
	//================================Start my SUS_MC simulation=================================
    // Weights WF = [(s = -0.5),(s=-0.4995),(s=-0.4990),(s=0.4985),(s=0.4980) ... , (s=0.3)] ---> say 0.8/0.0005 = 1600 components
    // Windows = [{(s = -0.5),(s = -0.4990)},{(s=-0.4995),(s=0.4985)}; {(s=-0.4990),(s=0.4980)};...;{(s=0.2990),(s=0.3)}] ---> 0.8/0.0005 = 1600 windows
    // !! the while loop start from a state that when S is in -0.5 !! NOT EMPTY!!!
	//===========================================================================================
	// std::tuple<int,int,int,int> Ad;
	// std::tuple<int,int,int,int> De;
    
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
		    size = nv+nh+nu;

    		aaccp = z*double(n0*n1*n2)/(double(nh+nv+nu+1.0)*double(length));
			daccp = (double(nh+nv+nu)*double(length))/(z*double(n0*n1*n2));

			prob = ((double) rand() / (RAND_MAX)); 

			// aaccp = (z*V)/((size+1.0)*K)*(exp(WF[int(w)] - WF[int(w-1)]));
			// daccp = (size*K)/(z*V)*(exp(WF[int(w-1)] - WF[int(w)]));	

			proba = min(1.0,aaccp);
			probd = min(1.0,daccp);

	        // ===========================Addition ===================================
			if(addordel == 0) 
			{
				AddSUS(s,prob,proba,w,DeltaS);
				// cout << (nu-(nv+nh)*0.5)/(nu+nv+nh)<<endl;
			}
			// ============================Deletion=============================
			else 
			{
				if (size != 0) // make sure there are rods to be del;
				{
					DelSUS(s,prob,probd,size,w,DeltaS);	
					// cout << std::get<0>(De) <<"  "<<std::get<1>(De)<<"  "<< std::get<2>(De)<<"  "<<std::get<3>(De)<<endl;
				}
			}	

			S = (nu-(nv+nh)*0.5)/(nu+nv+nh);

			// Now check S, only accept move when S is inside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
			if (-0.5+(w/2.0)*DeltaS <= S && S <= -0.5+(w+1.0)*DeltaS/2.0) // update the fu after each step.
			{
				fu++; // if at the upper window, update fu
				PH_w[w][size]++; // update the distribution/histogram of PH;
			}
			else if (-0.5+(w-1.0)*DeltaS/2.0 <= S && S <= -0.5+(w/2.0)*DeltaS)  // update the fl after each step.
			{
				fl++;//if at the lower window, update fl
				PL_w[w][size]++; // update the distribution/histogram of PL;
			}	
		}
		// cout << fu<< "  "<<fl<<"  "<<S<<endl;
		// =======================  if fu and fl != 0 Update the upper window side ================================
        if (fu!=0 && fl != 0 && (-0.5+(w/2.0)*DeltaS <= S))
        {

		    WF[w] = WF[w] + log(fu/fl);
	        // "linearly extrapolate" for WF[w+1] by using W[w] and WF[w-1]
	        WF[w+1] = WF[w];

	        cout << "# of Up rod: "<< nu<<"  # of Hor rod: "<<nh <<" # of Ver rod: " <<nv <<"  S = "<<S<<endl;
	        data_s << "# of Up rod: "<< nu<<"  # of Hor rod: "<<nh <<" # of Ver rod: " <<nv <<"  S = "<<S<<endl;

			// ======================= Print out the data into terminal =============================================		
			cout <<"Window: "<< w <<" : "<<"W("<<w<<" : lower) = "<< WF[w-1]<<" "<<"W("<<w<<" : Upper) = "<< WF[w] << endl;
			data_s <<"Window: "<< w <<" : "<<"W("<<w<<" : lower) = "<< WF[w-1]<<" "<<"W("<<w<<" : Upper) = "<< WF[w] << endl;

			// initial config determine the intial value of fu and fl

		    w++; // switch into the next window
        }
        // else reset PL_w[w] and PH_w[w], fu and fl and repeat the window simulation
        else
        {
        	for ( int i = 0;i<n0*n1*n2/length; i++ )
        	{
	        	PH_w[w][i] = 0;
	        	PL_w[w][i] = 0;
        	}
        }
	}

	double C=0.0; // the normalization constant for P_N


	// Calc the distribution P_w(N) by PH_w(N) and PL_w(N),Then Extrapolate the P_N by P_w and WF;
    double WFMAX = 0.0;
	for (int w = 1;w<N_window+1;w++)
	{
		if(WF[w] >WFMAX )
		{
			WFMAX = WF[w];
		}

		for(int i = 0; i<n0*n1*n2/length; i++)
		{
			P_w[w][i] = PL_w[w+1][i] + PH_w[w][i]; //Calc the distribution P_w(N)
			CN[w] += P_w[w][i];
		}
		//Normalize P_ws
		for (int i =0; i<n0*n1*n2/length; i++)
		{
			P_w[w][i] = P_w[w][i]/CN[w];
		}
	}

/*
    Extrapolate P_w_Z2 by assumming P_w_Z2[i] = P_w[i]*Z2^{i} and P_w(S) = C*exp(W(S))
    Then lnP_w_Z2[i] = lnP_w[i] + i *ln(Z2);
*/

 //    for (int w = 1;w<N_window+1;w++)
	// {
	//     vector<double> NORlnP_w, NORP_w; // an array of normalization constant for ln(P_w_Z2) and P_w_Z2
	//     NORP_w.resize(n0*n1*n2/length);
	//     NORlnP_w.resize(n0*n1*n2/length);

	// 	for(int i = 0; i<n0*n1*n2/length; i++)
	// 	{
	// 	    P_w_Z2[w][i] = log(P_w[i]) + i*log(9.33); // P_w_Z2 is actually ln(P_w_Z2) + C now;
	// 	    NORlnP_w[w] +=P_w_Z2[w][i];
	// 	}
	// 	//Normalize ln(P_w_Z2) + C 
	// 	for (int i =0; i<n0*n1*n2/length; i++)
	// 	{
	// 		P_w_Z2[w][i] = P_w_Z2[w][i]/NORlnP_w[w];
	// 		P_w_Z2[w][i] = exp(P_w_Z2[w][i]);
	// 		NORP_w[w] +=P_w_Z2[w][i];
	// 	}
	// 	for (int i =0; i<n0*n1*n2/length; i++)
	// 	{
	// 		P_w_Z2[w][i] = P_w_Z2[w][i]/NORP_w[w];
	// 	}
	// }

//Extrapolate the P_N by P_w and WF;

/*	
	    P(N) = C*Σ_{i=1}^{w_max} {P(i,N)e^{W(i)}} 
             = C*(P(1;N)*e^{W(1)} + P(2;N)*e^{W(2)} + ... +P(w_max;N)*e^{W(w_max)})
	P(N; Z2) = P(N)*Z2^{N} 
	         = C*Σ_{i=1}^{w_max} {P(i,N)e^{W(i)}}*Z2^{N}
  lnP(N; Z2) = lnC + ln(Σ_{i=1}^{w_max} {P(i,N)e^{W(i)}}) + ln(Z2^{N}) 
             = lnC + ln(P(1;N)*e^{W(1)} + P(2;N)*e^{W(2)} + ... +P(w_max;N)*e^{W(w_max)}) + N*ln(Z2)
             = lnC + ln(P(N)) + N*ln(Z2) -->but we can do this because P(N) has 0s and ln0 gives -inf

*/


	for(int k = 0; k<n0*n1*n2/length;k++) 
	{
		// Now calc the final P_N by using WF and P_w;
		for (int i = 1; i< N_window+1; i++)
		{
	        P_N[k] += P_w[i][k]*exp(WF[i]-WFMAX + 400);
	        // P_N_Z2[k] += P_w_Z2[i][k]*exp(WF[i]-1900);
		}
        C+=P_N[k];
    }


	for (int i = 1; i< N_window+1; i++)
	{
		sh<<WF[i]<<endl; // record WF into a text file;
	}

	double max = 0;
	double Coutsum = 0;

	for (int i =0; i<n0*n1*n2/length; i++)
	{
		P_N[i] = P_N[i]/C; // Normalized it
		ph<<P_N[i]<<endl;
	}

	// for(int i =0; i<n0*n1*n2/length; i++)
	// {
	// 	P_ext[i] = P_ext[i] - max + 500;
	// 	P_ext[i] = exp(P_ext[i]);
	// 	Coutsum += P_ext[i];
	// }

	// for(int i =0; i<n0*n1*n2/length; i++)
	// {
	// 	P_ext[i] = P_ext[i]/Coutsum;
	// 	// extrap <<P_ext[i]<<endl;
	// }


	ofstream myfile ("SUSWeight_function.txt");
	ofstream data ("Data.txt");
	// ofstream myfile2("P_w.txt");
	ofstream myfile3("P_N.txt");
	// ofstream myfile4("P_9_33.txt");

	// myfile2.precision(20);
	myfile3.precision(20);
	data.precision(20);
	myfile.precision(20);
	// myfile4.precision(20);

	string data4 = data_s.str();
	// string data3 = ph.str();
	string data5 = ph.str();
	string data2 = sh.str();
	// string data6 = extrap.str();

	data << data4;
	myfile << data2;
	// myfile2 << data3;
	myfile3 << data5;
	// myfile4 << data6;

	myfile.close();
	// myfile2.close();
	myfile3.close();
	data.close();
	// myfile4.close();



	return ;  
}










