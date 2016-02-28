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


std::tuple<int,int,int,int> MC::Add(Cells &s,double &prob,double &proba)
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
				    return std::make_tuple(x,y,z,o);					
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
				    return std::make_tuple(x,y,z,o);
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
				    return std::make_tuple(x,y,z,o);
				}
			}							
		}
    }
    return std::make_tuple(-1,-1,-1,-1);
}
void MC::Add_redo(Cells &s, std::tuple<int,int,int,int> Ad)
{
	int x,y,z,o;
	x = std::get<0>(Ad);
	y = std::get<1>(Ad);
	z = std::get<2>(Ad);
	o = std::get<3>(Ad);
	// delete the rod we just added;
	if(o == 0)
	{
		// --------------------- it's a vertical rod -----------------------
		for(int i = 0; i<length; i++)
		{
			// update the new config of cells
			s.getSquare(x,(y+i)%n1,z).setStatus(0);
		}
		// remove the target rod from the vector Rodlist;
		VRodlist.pop_back();
		av--;
		nv--;// cancel the earlier accumulate the # of ver rod;
	}
	
	else if (o == 1)
	{
		// --------------------- it's a horizontal rod -----------------------

		for(int i = 0; i<length; i++)
		{
			// update the new config of cells
			s.getSquare((x+i)%n0,y,z).setStatus(0);
		}
		// remove the target rod from the vector Rodlist;
		HRodlist.pop_back();		
		ah--;
		nh--;// cancel the earlier accumulate the # of hor rod;
	}
	else
	{
		// --------------------- it's a up rod -----------------------
		for(int i = 0; i<length; i++)
		{
			// update the new config of cells
			s.getSquare(x,y,(z+i)%n2).setStatus(0);
		}
		// remove the target rod from the vector Rodlist;
		URodlist.pop_back();
		au--;
		nu--;// cancel the earlier accumulate the # of hor rod;
	}	
}

std::tuple<int,int,int,int> MC::Del(Cells &s,double &prob,double &probd, double &size)
{		
	int x,y,z,indx;// the position of the target on the cells;
	x=y=z=indx=0;

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
			
			else if (indx < nv + nh && nh !=0)
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
			else if(nu !=0)
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
			return std::make_tuple(x,y,z,indx);					
		}
		return std::make_tuple(-1,-1,-1,-1); // if delesion failed, return -1,-1,-1,-1
	}
}

void MC::Del_redo(Cells &s,std::tuple<int,int,int,int> De)
{
	int x,y,z,indx,o;
	x = std::get<0>(De);
	y = std::get<1>(De);
	z = std::get<2>(De);
	indx = std::get<3>(De);

	if(indx < nv)
	{
		o = 0;
		HR rod(x,y,z,length,o);		
	//======================== Vertical ===============================
		// Do addition to put the deleted rod back;
		// push the rod back into the VRodlist;
		VRodlist.insert(VRodlist.begin() + indx,rod);
		// VRodlist.push_back(rod);
		nv++;
		dv--;
		// update new N, E and new config;
		for (int i = 0; i < length; i++)
		{	
			s.getSquare(x,(y+i)%n1,z).setStatus(1);
		}													
	}

	else if(indx < nv + nh)
	{
		o = 1;
		HR rod(x,y,z,length,o);		
    //======================= Horizontal  ============================
		// Do addition to put the deleted rod back;
		// push the rod into the HRodlist;
		HRodlist.insert(HRodlist.begin() + indx - nv,rod);
		// HRodlist.push_back(rod);
		nh++;
		dh--;
		// update new N, E and new config;
		for (int i = 0; i < length; i++)
		{
			s.getSquare((x+i)%n0,y,z).setStatus(1);
		}
	}
	else 
	{
		o = 2;
		HR rod(x,y,z,length,o);		
    //======================= Up  ============================
		// Do addition to put the deleted rod back;
		// push the rod into the URodlist;
		URodlist.insert(URodlist.begin()+indx-nv-nh,rod);
		// URodlist.push_back(rod);
		nu++;
		du--;
		// update new N, E and new config;
		for (int i = 0; i < length; i++)
		{
			s.getSquare(x,y,(z+i)%n2).setStatus(1);
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
* At the end, for Window i: we should have P_wi(N) = P_H_wi(N)+P_L_wi+1(N)
* and we obtained the P(N) by P(N) = C*sum_w1->wmax{P_wi(N)*exp{WF(i)}}
*
* According to R.L.C.Vink's paper, Delta S = 0.001 to 0.002 gives good result
*********************************************************************************
*/	

array<double,2000> MC::MCSUS()
{
	Cells s(n0,n1,n2,EMPTY,length);  //  setting the lattice;
	//================= declare instance variables ============================= 
	stringstream sh,ph,data_s;
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
    array<double,2000> WF; // NOTICED THAT THE # OF WEIGHTS IS NOT THE TOTAL NUMBER ANYMORE!!!!   --> 1D array
    // -----------------------The rest arrays are 2D arrays ------------------------------//
	vector<vector<double>> P_w; // an array that stores the histogram/distribution of ave # of N
	vector<vector<double>> PH_w; // an array that stores the histogram/distribution of Higher window in terms of # of N
	vector<vector<double>> PL_w; // an array that stores the histogram/distribution of Lower window in terms of # of N
	PH_w.resize(2000);
	PL_w.resize(2000);
	P_w.resize(2000);

	for(int i =0; i<2000; i++)
	{
		PH_w[i].resize(n0*n1*n2/length);
		PL_w[i].resize(n0*n1*n2/length);
		P_w[i].resize(n0*n1*n2/length);
		// cout<< "LOL: "<<i<<endl;
	}

    // ============= Do a GCMC to reach a initial state where S is in [0,deltaS] ==========
	while (S <= 0 || S>DeltaS )
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
	//================================Start my SUS_MC simulation=================================
    // Weights WF = [(s=0.0005),(s=0.001),(s=0.0015),(s=0.002) ... , (s=0.5)] ---> say 0.5/0.0005 = 1000 components
    // Windows = [{(s=0.0005),(s=0.001)}; {(s=0.001),(s=0.0015)};...;{(s=0.4995),(s=0.5)}] ---> 0.4995/0.0005 = 999 windows
    // !! the while loop start from a state that when S is in [0,deltaS] !! NOT EMPTY!!!
	//===========================================================================================
	std::tuple<int,int,int,int> Ad;
	std::tuple<int,int,int,int> De;
    
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

			prob = ((double) rand() / (RAND_MAX)); 

	        // ===========================Addition ===================================
			if(addordel == 0) 
			{
				Ad = Add(s,prob,proba);
				// cout << (nu-(nv+nh)*0.5)/(nu+nv+nh)<<endl;
			}
			// ============================Deletion=============================
			else 
			{
				if (size != 0) // make sure there are rods to be del;
				{
					De = Del(s,prob,probd,size);	
					// cout << std::get<0>(De) <<"  "<<std::get<1>(De)<<"  "<< std::get<2>(De)<<"  "<<std::get<3>(De)<<endl;
				}
			}	

			S = (nu-(nv+nh)*0.5)/(nu+nv+nh);

			// Now check S, only accept move when S is inside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
			// if ((w-1.0)*DeltaS/2.0<= S && S <= (w+1.0)*DeltaS/2.0)
			// {
			if ((w/2.0)*DeltaS <= S && S <= (w+1.0)*DeltaS/2.0) // update the fu after each step.
			{
				fu++; // if at the upper window, update fu
				PH_w[w][size]++; // update the distribution/histogram of PH;
			}
			else if ((w-1.0)*DeltaS/2.0 <= S && S < (w/2.0)*DeltaS)  // update the fl after each step.
			{
				fl++;//if at the lower window, update fl
				PL_w[w][size]++; // update the distribution/histogram of PL;
			}	
			// }
			else
			{	
				// the case that S is outside [(w-1.0)*DeltaS/2.0 , (w+1.0)*DeltaS/2.0)]
				if(addordel == 0) 
				{
					if(std::get<0>(Ad) !=-1 && std::get<1>(Ad) !=-1 && std::get<2>(Ad) !=-1 && std::get<3>(Ad) !=-1)
					{ 
						//make sure we did Add successfully before;
						Add_redo(s,Ad);
					}
				}
				else 
				{
					if(std::get<0>(De) !=-1 && std::get<1>(De) !=-1 && std::get<2>(De) !=-1 && std::get<3>(De) !=-1)
					{
						//make sure we did Del successfully before;
						Del_redo(s,De);
					}
				}	
			} 
		}
		// =======================  if fu and fl != 0 Update the upper window side ================================
        if (fu!=0 && fl != 0)
        {
		    WF[w] = WF[w] + log(fu/fl);
	        // "linearly extrapolate" for WF[w+1] by using W[w] and WF[w-1]
	        WF[w+1] = WF[w];

			// Store the distribution of PH_w(N) and PL_w(N) and update P_w(N)
			for(int i =0;i<n0*n1*n2/length; i++)
			{
				P_w[w][i]= PL_w[w][i] + PH_w[w-1][i];
			}
	        cout << "# of Up rod: "<< nu<<"  # of Hor rod: "<<nh <<" # of Ver rod: " <<nv <<"  S = "<<S<<endl;
	        data_s << "# of Up rod: "<< nu<<"  # of Hor rod: "<<nh <<" # of Ver rod: " <<nv <<"  S = "<<S<<endl;

			// ======================= Print out the data into terminal =============================================		
			cout <<"Window: "<< w <<" : "<<"W("<<w<<" : lower) = "<< WF[w-1]<<" "<<"W("<<w<<" : Upper) = "<< WF[w] << endl;
			data_s <<"Window: "<< w <<" : "<<"W("<<w<<" : lower) = "<< WF[w-1]<<" "<<"W("<<w<<" : Upper) = "<< WF[w] << endl;

			// initial config determine the intial value of fu and fl

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

	for(int i = 1; i< N_window+1; i++) 
	{
		for(int k = 0; k<n0*n1*n2/length; k++)
		{
			ph<< P_w[i][k] <<" ";
		}
		ph<<endl;

		sh<<WF[i]<<endl;
	}

	ofstream myfile ("SUSWeight_function.txt");
	ofstream data ("Data.txt");
	ofstream myfile2("P_w.txt");
	myfile2.precision(20);
	data.precision(20);
	myfile.precision(20);
	string data4 = data_s.str();
	string data3 = ph.str();
	string data2 = sh.str();
	data << data4;
	myfile << data2;
	myfile2 << data3;
	myfile.close();
	myfile2.close();
	data.close();

	return WF;  
}





