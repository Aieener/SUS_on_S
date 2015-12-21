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


// void MC::Add(Cells &s,double &prob,double &proba,double &V,double &DeltaS,double &w) //Revised for SUS on S!!!
// {
// 	int x,y,z,o; // pick a random position and orientation for the HR to be added;
// 	x = rand()%n0; // x range[0,n0-1]
// 	y = rand()%n1; // y range[0,n1-1]
// 	z = rand()%n2; // z range[0,n2-1]
// 	o = rand()%3  ; // 0 range {0,1,2}
// 	double S_ = 0.0;

// 	if(s.getSquare(x,y,z).isEmpty()) // if it's open, try to do Addition;
// 	{
// 		HR rod(x,y,z,length,o);

// 		if(o == 0)
// 		{
// 		//======================== Vertical ===============================
// 			// the vertical case
// 			int counter = 0;

// 			for (int j = 0; j < length-1; j++)
// 			{
// 				// check if the vertical space is wide open
// 				if(s.getSquare(x,(y+j+1)%n1,z).isOccupied())
// 				{
// 					counter++;
// 					break;
// 				}						
// 			}
// 			if (counter == 0)
// 			{
// 				if(prob<=proba)
// 				{
// 					av++;
// 					nv++;// accumulate the # of ver rod;
// 				    S_ = abs((nu/V-(nv/V+nh/V)/2)/(nv/V+nh/V+nu/V));
// 				    cout << S_<< endl;	
// 					if(S_>= (w-1.0)*DeltaS/2.0 && S_<(w+1.0)*DeltaS/2.0)// CHECK WINDOW BOUND
// 					{					    
// 						// Do addition;
// 						// push the new rod into the VRodlist;
// 						VRodlist.push_back(rod);
// 						// update new N, E and new config;
// 						for (int i = 0; i < length; i++)
// 						{	
// 							s.getSquare(x,(y+i)%n1,z).setStatus(1);
// 						}
// 					}	

// 					else
// 					{
// 						// Addition out of bound of the window, Do nothing and take the av and nv back to their old value
// 						av--;
// 						nv--; 
// 					}				
// 				}		
// 			}									
// 		}

// 		else if(o == 1)
// 		{
//         //======================= Horizontal  ============================
// 			int counter = 0;
// 			for (int j = 0; j< length-1 ; j++)
// 			{
// 				// check if the horizontal space is wide open
// 				if(s.getSquare((x+1+j)%n0,y,z).isOccupied())
// 				{
// 					counter++;
// 					break;
// 				}							
// 			}
// 			if (counter == 0)
// 			{
// 				if(prob<= proba)
// 				{
// 					ah++;
// 					nh++;// accumulate the # of hor rod;
// 					S_ = abs((nu/V-(nv/V+nh/V)/2)/(nv/V+nh/V+nu/V));	
// 					if(S_>= (w-1.0)*DeltaS/2.0 && S_<(w+1.0)*DeltaS/2.0)// CHECK WINDOW BOUND
// 					{					 
// 						//Do addition;
// 						//push the new rod into the HRodlist;
// 						HRodlist.push_back(rod);
// 						// update new N, E and new config;
// 						for (int i = 0; i < length; i++)
// 						{
// 							s.getSquare((x+i)%n0,y,z).setStatus(1);
// 						}
// 					}
// 					else
// 					{
// 						// Addition out of bound of the window, Do nothing and take the av and nv back to their old value
// 						ah--;
// 						nh--; 
// 					}	
// 				}
// 			}
// 		}

// 		else 
// 		{
//         //======================= Up  ============================
// 			int counter = 0;
// 			for (int j = 0; j< length-1 ; j++)
// 			{
// 				// check if the horizontal space is wide open
// 				if(s.getSquare(x,y,(z+j+1)%n2).isOccupied())
// 				{
// 					counter++;
// 					break;
// 				}							
// 			}
// 			if (counter == 0)
// 			{
// 				if(prob<= proba)
// 				{
// 					au++;
// 					nu++;// accumulate the # of hor rod;
// 					S_ = abs((nu/V-(nv/V+nh/V)/2)/(nv/V+nh/V+nu/V));	
// 					if(S_>= (w-1.0)*DeltaS/2.0 && S_<(w+1.0)*DeltaS/2.0)// CHECK WINDOW BOUND
// 					{			
// 						//Do addition;
// 						//push the new rod into the HRodlist;
// 						URodlist.push_back(rod);
// 						// update new N, E and new config;
// 						for (int i = 0; i < length; i++)
// 						{
// 							s.getSquare(x,y,(z+i)%n2).setStatus(1);
// 						}
// 					}
// 					else
// 					{
// 						// Addition out of bound of the window, Do nothing and take the av and nv back to their old value
// 						au--;
// 						nu--; 
// 					}	
// 				}
// 			}							
// 		}
//     }
// }

// void MC::Del(Cells &s,double &prob,double &probd, double &size,double &V,double &DeltaS, double &w) //Revised for SUS on S!!!
// {
// 	double S_ = 0.0;
// 	if(nv + nh + nu >0)// make sure there are Vertical rod;
// 	{
// 		int indx; // pick a random index from the Rodlist;
// 		indx = rand()%int(nv+nh+nu);

// 		//remove Rodlist[indx];
// 		int x,y,z;// the position of the target on the cells;

// 		if(prob <= probd)
// 		{					
// 			if(indx < nv)
// 			{
// 				x = VRodlist[indx].getX();
// 				y = VRodlist[indx].getY();
// 				z = VRodlist[indx].getZ();	
// 				// --------------------- it's a vertical rod -----------------------
// 				nv--;// substract the # of ver rod;
// 				dv++;
// 				S_ = abs((nu/V-(nv/V+nh/V)/2)/(nv/V+nh/V+nu/V));

// 				if(S_>= (w-1.0)*DeltaS/2.0 && S_<(w+1.0)*DeltaS/2.0)// CHECK WINDOW BOUND
// 				{
// 					// Do Deletion
// 					for(int i = 0; i<VRodlist[indx].getLength(); i++)
// 					{
// 						// update the new config of cells
// 						s.getSquare(x,(y+i)%n1,z).setStatus(0);
// 					}
// 					// remove the target rod from the vector Rodlist;
// 					VRodlist.erase(VRodlist.begin() + indx);
// 				}
// 				else
// 				{
// 					//Deletion out of bound; recover nv and dv and do nothing
// 					nv++;
// 					dv--;
// 				}
// 			}
			
// 			else if (indx < nv + nh)
// 			{
// 				nh--;// substract the # of ver rod;
// 				dh++;
// 				S_ = abs((nu/V-(nv/V+nh/V)/2)/(nv/V+nh/V+nu/V));
// 				x = HRodlist[indx - nv].getX();
// 				y = HRodlist[indx - nv].getY();
// 				z = HRodlist[indx - nv].getZ();	
// 				// --------------------- it's a horizontal rod -----------------------
// 				if(S_>= (w-1.0)*DeltaS/2.0 && S_<(w+1.0)*DeltaS/2.0)// CHECK WINDOW BOUND
// 				{
// 					// Do Deletion
// 					for(int i = 0; i<HRodlist[indx-nv].getLength(); i++)
// 					{
// 						// update the new config of cells
// 						s.getSquare((x+i)%n0,y,z).setStatus(0);
// 					}
// 					// remove the target rod from the vector Rodlist;
// 					HRodlist.erase(HRodlist.begin() + indx - nv);
// 				}
// 				else
// 				{
// 					//Deletion out of bound; recover nv and dv and do nothing
// 					nh++;
// 					dh--;
// 				}
// 			}
// 			else
// 			{
// 				nu--;// substract the # of ver rod;
// 				du++;
// 				S_ = abs((nu/V-(nv/V+nh/V)/2)/(nv/V+nh/V+nu/V));
// 				if(S_>= (w-1.0)*DeltaS/2.0 && S_<(w+1.0)*DeltaS/2.0)// CHECK WINDOW BOUND
// 				{
// 					// Do Deletion
// 					x = URodlist[indx - nv - nh].getX();
// 					y = URodlist[indx - nv - nh].getY();
// 					z = URodlist[indx - nv - nh].getZ();	
// 					// --------------------- it's a up rod -----------------------
// 					for(int i = 0; i<URodlist[indx-nv-nh].getLength(); i++)
// 					{
// 						// update the new config of cells
// 						s.getSquare(x,y,(z+i)%n2).setStatus(0);
// 					}
// 					// remove the target rod from the vector Rodlist;
// 					URodlist.erase(URodlist.begin()+indx-nv-nh);
// 				}
// 				else
// 				{
// 					//Deletion out of bound; recover nv and dv and do nothing
// 					nu++;
// 					du--;
// 				}
// 			}
// 		}
// 	}	
// }

void MC::Add(Cells &s,array<double,100000> &WF,double &prob,double &size,double &S,double &deltaS)
{
	double aaccp,proba;      // the acceptance probabilities: 
	double V = double(n0*n1*n2);    // the total lattice size
	double K = double(length); //
	int x,y,z,o; // pick a random position and orientation for the HR to be added;
	x = rand()%n0; // x range[0,n0-1]
	y = rand()%n1; // y range[0,n1-1]
	z = rand()%n2; // z range[0,n2-1]
	o = rand()%3  ; // 0 range {0,1,2}
    int j_0 = S/(deltaS/2); // the initial state index
    // cout << S<<"  "<<j_0<<endl;

    int j_1; // the state after movement

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
				av++;
				nv++;// accumulate the # of ver rod;
				j_1 = S/(deltaS/2);
				aaccp = (z*V)/((size+1.0)*K)*(exp(WF[j_1] - WF[j_0])); // Have to revise!
				proba = min(1.0,aaccp);

				if(prob<=proba)
				{
					// Do addition;
					// push the new rod into the VRodlist;
					VRodlist.push_back(rod);
					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{	
						s.getSquare(x,(y+i)%n1,z).setStatus(1);
					}					
				}
				else
				{ // movement failed, refine av and nv
					av--;
					nv--;// accumulate the # of ver rod;
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
				ah++;
				nh++;// accumulate the # of hor rod;
				j_1 = S/(deltaS/2);
				aaccp = (z*V)/((size+1.0)*K)*(exp(WF[j_1] - WF[j_0])); // Have to revise!
				proba = min(1.0,aaccp);

				if(prob<= proba)
				{
					//Do addition;
					//push the new rod into the HRodlist;
					HRodlist.push_back(rod);
					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{
						s.getSquare((x+i)%n0,y,z).setStatus(1);
					}
				}
				else
				{ // movement failed, refine av and nv
					ah--;
					nh--;// accumulate the # of ver rod;
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
				au++;
				nu++;// accumulate the # of hor rod;
				j_1 = S/(deltaS/2);
				aaccp = (z*V)/((size+1.0)*K)*(exp(WF[j_1] - WF[j_0])); // Have to revise!
				proba = min(1.0,aaccp);
				if(prob<= proba)
				{
					//Do addition;
					//push the new rod into the HRodlist;
					URodlist.push_back(rod);

					// update new N, E and new config;
					for (int i = 0; i < length; i++)
					{
						s.getSquare(x,y,(z+i)%n2).setStatus(1);
					}
				}
				else
				{ // movement failed, refine av and nv
					au--;
					nu--;// accumulate the # of ver rod;
				}	
			}							
		}
    }
}

void MC::Del(Cells &s,array<double,100000> &WF,double &prob, double &size,double &S,double &deltaS)
{
	double daccp,probd;      // the acceptance probabilities:
	double V = double(n0*n1*n2);    // the total lattice size
	double K = double(length); // 
    int j_0 = S/(deltaS/2); // the initial state index
    int j_1; // the state after movement

	if(nv + nh + nu >0)// make sure there are Vertical rod;
	{
		int indx; // pick a random index from the Rodlist;
		indx = rand()%int(nv+nh+nu);

		//remove Rodlist[indx];
		int x,y,z;// the position of the target on the cells;

		if(indx < nv)
		{
			x = VRodlist[indx].getX();
			y = VRodlist[indx].getY();
			z = VRodlist[indx].getZ();	
			// --------------------- it's a vertical rod -----------------------
			nv--;// substract the # of ver rod;
			dv++;
			j_1 = S/(deltaS/2);			
			daccp = (size*K)/(z*V)*(exp(WF[j_1] - WF[j_0]));	//Have to revise!
			probd = min(1.0,daccp);

			if(prob <= probd)
			{
				//Do deletion
				for(int i = 0; i<VRodlist[indx].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare(x,(y+i)%n1,z).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				VRodlist.erase(VRodlist.begin() + indx);
			}
			else
			{
				// deletion failed, refine nv and dv
				nv++;
				dv--;
			}
		}
		
		else if (indx < nv + nh)
		{
			x = HRodlist[indx - nv].getX();
			y = HRodlist[indx - nv].getY();
			z = HRodlist[indx - nv].getZ();	
			// --------------------- it's a horizontal rod -----------------------
			nh--;// substract the # of ver rod;
			dh++;
			j_1 = S/(deltaS/2);			
			daccp = (size*K)/(z*V)*(exp(WF[j_1] - WF[j_0]));	//Have to revise!
			probd = min(1.0,daccp);
			if(prob <= probd)
			{
				for(int i = 0; i<HRodlist[indx-nv].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare((x+i)%n0,y,z).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				HRodlist.erase(HRodlist.begin() + indx - nv);				
			}
			else
			{
				nh++;
				dh--;
			}
		}
		else
		{
			x = URodlist[indx - nv - nh].getX();
			y = URodlist[indx - nv - nh].getY();
			z = URodlist[indx - nv - nh].getZ();	
			// --------------------- it's a up rod -----------------------
			nu--;// substract the # of ver rod;
			du++;
			j_1 = S/(deltaS/2);			
			daccp = (size*K)/(z*V)*(exp(WF[j_1] - WF[j_0]));	//Have to revise!
			probd = min(1.0,daccp);
			if(prob <= probd)
			{
				for(int i = 0; i<URodlist[indx-nv-nh].getLength(); i++)
				{
					// update the new config of cells
					s.getSquare(x,y,(z+i)%n2).setStatus(0);
				}
				// remove the target rod from the vector Rodlist;
				URodlist.erase(URodlist.begin()+indx-nv-nh);				
			}
			else
			{
				nu++;
				du--;
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
	double prob;               // the prob to decide either accept add/del;
	double V = double(n0*n1*n2);    // the total lattice size

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


	//================================Start my SUS_MC simulation=================================
    //Weights WF = [(s=0.0005),(s=0.001),(s=0.0015),(s=0.002) ... , (s=0.5)] ---> say 0.5/0.0005 = 1000 components
    //Windows = [{(s=0.0005),(s=0.001)}; {(s=0.001),(s=0.0015)};...;{(s=0.4995),(s=0.5)}] ---> 0.4995/0.0005 = 999 windows
	//===========================================================================================
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
		    double S = 0;
		    if(nu+nv+nh != 0)
		    {
		    	S = (nu-(nv+nh)*0.5)/(nu+(nv+nh)*0.5);
		    } 
			prob = ((double) rand() / (RAND_MAX)); 

	        // ===========================Addition ===================================
			if(addordel == 0) 
			{
				//For each window, N is freely fluctuated		
				//try to Do Addition;
				Add(s,WF,prob,size,S,DeltaS);
			}
			// ============================Deletion=============================
			else 
			{
				//For each window, N is freely fluctuated
				//Do deletion;
				Del(s,WF,prob,size,S,DeltaS);	
			}	

			if ((w/2.0)*DeltaS <= S && S < (w+1.0)*DeltaS/2.0) // update the fu after each step.
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
		}
		// =======================  if fu and fl != 0 Update the upper window side ================================
        if (fu!=0 && fl != 0)
        {
		    WF[w] = WF[w] + log(fu/fl);
	        // "linearly extrapolate" for WF[w+1] by using W[w] and WF[w-1]
	        // WF[w+1] = WF[w];
	        WF[w+1] = 2*WF[w] - WF[w-1];

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





