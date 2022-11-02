#include<iostream>
#include<istream>
#include<stdlib.h>
#include<string>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<cmath>
#include<vector>
#include<algorithm>
using namespace std;

// Code adapted from Emran Reshid's work.

//Power function
//----------------------------------------------------------------------------------------------------------------------------------------------//

double power(double val, int exp)	{
	double x = val*val;
	return x;
}

//Checking that system traps water
//----------------------------------------------------------------------------------------------------------------------------------------------//
bool waterTrappingCheck(double xcoordinate[], double ycoordinate[], double zcoordinate[], string ResidueID[], string MoleculeID[], double IdentifierNumber[], int numInputs, double dimensions[])	{
	//Binning code adapted from Kassandra Schmidt's work.
	int* beadNumber = new int[numInputs];

	for (int i = 0; i < numInputs; i++)	{
		beadNumber[i] = (int)IdentifierNumber[i];
	}

	//Ensuring box is big enough to have more than one bin
	double x_distance2 = dimensions[0] * dimensions[0];
    double y_distance2 = dimensions[1] * dimensions[1];
    double z_distance2 = dimensions[2] * dimensions[2];

    double box_diagonal_distance = sqrt(power(x_distance2,2) + power(y_distance2,2) + power(z_distance2,2));

    //Cut off distance, 1.2 nm diagonal
    double r_max = 1.2;

    //If box is not big enough, every particle is a neighbor to every other particle
    if (box_diagonal_distance < r_max)	{
    	return false;
    }

    //Otherwise, determine which particles are interacting with which other particles
    else	{

    	//Calculate number of bins. Must be whole number (rounds down).
    	int number_of_x_bins = int(dimensions[0] / r_max);
    	int number_of_y_bins = int(dimensions[1] / r_max);
    	int number_of_z_bins = int(dimensions[2] / r_max);

    	//Initialize bin list array
    	vector<int> binList [number_of_x_bins][number_of_y_bins][number_of_z_bins];
    	//vector<vector<vector<vector<int>>>> binList(number_of_x_bins, vector<vector<int>>(number_of_y_bins, vector<int>(number_of_z_bins)));
    	//vector<vector<vector<vector<int> > > > binList(3, vector<vector<vector<int> > >(4, vector<vector<int> >(5)));


    	//Adjusts bin widths to ensure all particles are in a bin.
    	double x_bin_width = dimensions[0] / number_of_x_bins;
    	double y_bin_width = dimensions[1] / number_of_y_bins;
    	double z_bin_width = dimensions[2] / number_of_z_bins;

    	for (int i = 0; i < numInputs; i++)	{
    		//cout<<i<<endl;
    		//Determining the right bins for x, y, z coordinates, accounting for particles on edges of bins and the box.
    		int particle_x_bin_number = int(xcoordinate[i] / x_bin_width);
    		if (particle_x_bin_number >= number_of_x_bins)	{
    			particle_x_bin_number = number_of_x_bins - 1;
    		}
    		//cout<<particle_x_bin_number<<endl;
    		//cout<<xcoordinate[i]<<endl;

    		int particle_y_bin_number = int(ycoordinate[i] / y_bin_width);
    		if (particle_y_bin_number >= number_of_y_bins)	{
    			particle_y_bin_number = number_of_y_bins - 1;
    		}
    		//cout<<particle_y_bin_number<<endl;
    		//cout<<ycoordinate[i]<<endl;

    		int particle_z_bin_number = int(zcoordinate[i] / z_bin_width);
    		if (particle_z_bin_number >= number_of_z_bins)	{
    			particle_z_bin_number = number_of_z_bins - 1;
    		}
    		//cout<<particle_z_bin_number<<endl;
    		//cout<<zcoordinate[i]<<endl;

    		//Adding the particle's atom ID to its corresponding bin.
    		binList[particle_x_bin_number][particle_y_bin_number][particle_z_bin_number].push_back(i); //Push the index of the bead into the bin list

    	}

    	//Checking that they are placed in the right position
    	for (int i = 0; i < number_of_x_bins; i++)	{
    		for (int j = 0; j < number_of_y_bins; j++)	{
    			for (int k = 0; k < number_of_z_bins; k++)	{
    				for (int l = 0; l < binList[i][j][k].size(); l++)	{
    					//cout << "x: " << i << " y: " << j << " z: " << k << endl;
    					//cout << binList[i][j][k][l] <<endl;
    				}
    			}
    		}
    	}

    	//Trajectory position coordinates in an array (updated as we proceed along a vector)
    	double trajectoryPosition[3] = {dimensions[0]/2,dimensions[1]/2,dimensions[2]/2}; //Initialized to the center of the box
    	double increment = 0.3; //Size to increment across the trajectory
    	double checkingRadius = 1.2; //Radius at which we scan for particles at a particle point along the vector
    	//int count = 0; //Count for the number of particles checked

    	//Proceed across a trajectory, check for water trapping 
    	for (int x = 0; x<3; x++)	{

    		vector<int>particleTypes; //List containing all of the particle types (peptide or water)

			bool inBox = true; //Boolean flag that indicates whether the point along the trajectory is still in the box
			trajectoryPosition[0] = dimensions[0]/2;
			trajectoryPosition[1] = dimensions[1]/2;
			trajectoryPosition[2] = dimensions[2]/2;
			cout<<dimensions[0]<<dimensions[1]<<dimensions[2]<<endl;

	    	while (inBox)	{
	    		int closestParticle = -1; //The particle closest to the point along the vector
	    		double closestParticleDistance = 1000000;

	    		//Look at the 27 relevant bins 
	    		for (int i = -1; i < 2; i++)	{
	    			for (int j = -1; j < 2; j++)	{
	    				for (int k = -1; k < 2; k++)	{
	    					//Find current bin
	    					int particle_x_bin_number = int(trajectoryPosition[0] / x_bin_width);
				    		if (particle_x_bin_number >= number_of_x_bins)	{
				    			particle_x_bin_number = number_of_x_bins - 1;
				    		}
				    		//cout << particle_x_bin_number << endl;

				    		int particle_y_bin_number = int(trajectoryPosition[1] / y_bin_width);
				    		if (particle_y_bin_number >= number_of_y_bins)	{
				    			particle_y_bin_number = number_of_y_bins - 1;
				    		}
				    		//cout << particle_y_bin_number << endl;

				    		int particle_z_bin_number = int(trajectoryPosition[2] / z_bin_width);
				    		if (particle_z_bin_number >= number_of_z_bins)	{
				    			particle_z_bin_number = number_of_z_bins - 1;
				    		}
				    		//cout << particle_z_bin_number << endl;

				    		//Adjusting bin number (so that we eventually look at all 27 bins)
				    		particle_x_bin_number = particle_x_bin_number + i;
				    		particle_y_bin_number = particle_y_bin_number + j;
				    		particle_z_bin_number = particle_z_bin_number + k;

				    		//Periodic boundary conditions
				    		if (particle_x_bin_number >= number_of_x_bins)	{
				    			particle_x_bin_number = 0;
				    		}
				    		else if (particle_x_bin_number < 0)	{
				    			particle_x_bin_number = number_of_x_bins -1;
				    		}

				    		if (particle_y_bin_number >= number_of_y_bins)	{
				    			particle_y_bin_number = 0;
				    		}
				    		else if (particle_y_bin_number < 0)	{
				    			particle_y_bin_number = number_of_y_bins -1;
				    		}

				    		if (particle_z_bin_number >= number_of_z_bins)	{
				    			particle_z_bin_number = 0;
				    		}
				    		else if (particle_z_bin_number < 0)	{
				    			particle_z_bin_number = number_of_z_bins -1;
				    		}

				    		for (int l = 0; l < binList[particle_x_bin_number][particle_y_bin_number][particle_z_bin_number].size(); l++)	{
				    			//cout<<xcoordinate[binList[particle_x_bin_number][particle_y_bin_number][particle_z_bin_number][l]]<<endl;
				    			double distance = sqrt(power((xcoordinate[binList[particle_x_bin_number][particle_y_bin_number][particle_z_bin_number][l]]) - trajectoryPosition[0],2) + power((ycoordinate[binList[particle_x_bin_number][particle_y_bin_number][particle_z_bin_number][l]]) - trajectoryPosition[1],2) + power((zcoordinate[binList[particle_x_bin_number][particle_y_bin_number][particle_z_bin_number][l]-1]) - trajectoryPosition[2],2));
				    			//cout<<"dist"<<distance<<endl;
				    			if ((distance < checkingRadius) && (distance < closestParticleDistance))	{
				    				closestParticle = binList[particle_x_bin_number][particle_y_bin_number][particle_z_bin_number][l];
				    				closestParticleDistance = distance;
				    			}
		    				}
				    	}
				    }
	    		}

	    		//Print statements to visualize the water trapping iterations
	    		//cout << IdentifierNumber[closestParticle] << endl;
	    		if (MoleculeID[closestParticle] == "S")	{
	    			cout << "SOLVENT" << endl;
					particleTypes.push_back(0);
				}
				else	{
					cout << "PEPTIDE" << endl;
					particleTypes.push_back(1);
				}

				//Iterating in the appropriate directions (x, y, z directions, respectively)
				if (x==0)	{
					double temp = trajectoryPosition[0];
					trajectoryPosition[0] = temp + increment;
				}
				else if (x==1)	{
					double temp = trajectoryPosition[1];
					trajectoryPosition[1] = temp + increment;
				}
				else	{
					double temp = trajectoryPosition[2];
					trajectoryPosition[2] = temp + increment;
				}

				if (abs(trajectoryPosition[0]) >= dimensions[0] || abs(trajectoryPosition[1]) >= dimensions[1] || abs(trajectoryPosition[2]) >= dimensions[2])	{
					//cout << trajectoryPosition[0] << endl;
					inBox = false;
				}
	    	}

	    	int numParticleChanges = 0; //Number of times it changes particle types

	    	/*for (int i = 0; i < particleTypes.size(); i++)	{
	    		cout << particleTypes[i] << endl;
	    	}*/

	    	//Going through the particle types list, ensuring that it starts with water, changes to peptide, then changes to water
	    	int numConsecutive = 0; //Number of consecutive same particle types
	    	bool firstParticle = true; //Boolean to use to check that the first particle is solvent
	    	int i = 0;

	    	while (i < particleTypes.size())	{
	    		//cout<<i<<endl;
	    		numConsecutive++;
	    		if (i > 0)	{
	    			if (particleTypes[i] != particleTypes[i-1])	{
	    				//cout << i << endl;

	    				//To account for any possible peptide particles that are randomly scattered
		    			if (numConsecutive > 3 && particleTypes[i+1] == particleTypes[i] && particleTypes[i+2] == particleTypes[i])	{ 
		    				//cout<<"hello"<<endl;
		    				//cout<<particleTypes[i]<<endl;
		    				//cout<<particleTypes[i-1]<<endl;
		    				numParticleChanges++;
		    			}
		    			numConsecutive = 0;
		    		}
	    		}

	    		if (firstParticle)	{
	    			if ((numConsecutive > 1) && (firstParticle))	{
		    			if (particleTypes[i] != 0)	{
		    				cout << i << endl;
				    		cout << "First particle is not solvent." << endl;
				    		return false;
				    	}
				    	firstParticle = false;
		    		}
	    		}
	    		i=i+1;
	    	}

	    	cout << "NUM PARTICLE CHANGES: " << numParticleChanges << endl;
	    	if (numParticleChanges == 2 || numParticleChanges == 4)	{
	    		return true;
	    	}
	    	std::vector<int>().swap(particleTypes);
	    }
	    delete[]beadNumber;

	    for (int i = 0; i < number_of_x_bins; i++)	{
	    	for (int j = 0; j < number_of_y_bins; j++)	{
	    		for (int k = 0; k < number_of_z_bins; k++)	{
	    			std::vector<int>().swap(binList[i][j][k]);
	    		}
	    	}
	    }

    	return false;
    }
}

//Nanovesicle characterization
//----------------------------------------------------------------------------------------------------------------------------------------------//
bool Nanovesicle(double xcoordinate[], double ycoordinate[], double zcoordinate[], string ResidueID[], string MoleculeID[], double IdentifierNumber[], int numInputs, double dimensions[])	{

	double a1;//1 represents (a,0,0)
	double findingacombo;
	double a1combo;
	int a1comboplaceholder;

	double avgDistYZ = 0;

	for(int r=0;r<numInputs;r++)//finds closest point to (a,0,0)
	{
		//cout<< ycoordinate[r] << " " << zcoordinate[r] << endl;
		findingacombo=sqrtl(power(fabs(ycoordinate[r]),2)+power(fabs(zcoordinate[r]),2));//finds the least positive distance to the center of mass for y and z coordinates
		/*if (r==0)
		{
			a1combo=findingacombo;
			a1comboplaceholder=r;
		}
		    
		if (a1combo>findingacombo)
		{
			a1combo=findingacombo;
			a1comboplaceholder=r;
		}*/
		avgDistYZ+=findingacombo;
		//cout<<"avgDistYZ: " << avgDistYZ << endl;
	}
	cout<<"avgDistYZ: " << avgDistYZ << endl << numInputs << endl;;
	avgDistYZ = avgDistYZ/numInputs;

	//a1=xcoordinate[a1comboplaceholder];

	double b1;//1 represents (0,b,0)
	double findingbcombo;
	double b1combo;
	int b1comboplaceholder;

	double avgDistXZ = 0;

	for(int r=0;r<numInputs;r++)//finds closest point to (0,b,0)
	{
		findingbcombo=sqrtl(power(fabs(xcoordinate[r]),2)+power(fabs(zcoordinate[r]),2));//finds the least positive distance to the center of mass for x and z coordinates
		/*if (r==0)
		{
			b1combo=findingbcombo;
			b1comboplaceholder=r;
		}
		    
		if (b1combo>findingbcombo)
		{
			b1combo=findingbcombo;
			b1comboplaceholder=r;
		}*/
		avgDistXZ+=findingbcombo;
		//cout<<"avgDistXZ: " << avgDistXZ<< endl;
	}

	avgDistXZ = avgDistXZ/numInputs;

	//b1=ycoordinate[b1comboplaceholder];

	double c1;//1 represents (0,0,c)
	double findingccombo;
	double c1combo;
	int c1comboplaceholder;
	double avgDistXY = 0;

	for(int r=0;r<numInputs;r++)//finds closest point to (0,0,c)
	{
		findingccombo=sqrtl(power(fabs(xcoordinate[r]),2)+power(fabs(ycoordinate[r]),2));//finds the least positive distance to the center of mass for x and y coordinates
		/*if (r==0)
		{
			c1combo=findingccombo;
			c1comboplaceholder=r;
		}
		    
		if (c1combo>findingccombo)
		{
			c1combo=findingccombo;
			c1comboplaceholder=r;
		}*/
		avgDistXY+=findingccombo;
		//cout<<"avgDistXY: " << avgDistXY << endl;
	}
	avgDistXY = avgDistXY/numInputs;

	//c1=zcoordinate[c1comboplaceholder];

	double VesicleSuccessCounter=0;
	double equationofvesicle;
	for(int r=0;r<numInputs;r++)
	{
		equationofvesicle=(power(xcoordinate[r],2)/power(avgDistYZ,2))+(power(ycoordinate[r],2)/power(avgDistXZ,2))+(power(zcoordinate[r],2)/power(avgDistXY,2));
		//cout<<"Equation of vesicle: "<< equationofvesicle<< endl;
		if ((0.5<equationofvesicle)&&(equationofvesicle<2))//acceptable range can be changed; currently : .8<=(x^2/a^2)+(y^2/b^2)+(z^2/c^2)<=1.2
		{
			VesicleSuccessCounter+=1;
		}
	}

	cout<<"Avg YZ dist: " << avgDistYZ << endl;
	cout<<"Avg XZ dist: " << avgDistXZ << endl;
	cout<<"Avg XY dist: " << avgDistXY << endl;
	
	cout<<"Percentage of Points that fits Nanovesicle structure : "<<VesicleSuccessCounter/(numInputs)*100<<" %"<<endl;

	cout<<endl;
	if (VesicleSuccessCounter/(numInputs)*100 > 70)	{ //Set 70% as the cut off mark
		return true;
	}
	return false;
}

int main()
{
	string FileName="peptideTextFile.txt";//Insert File name here; Reason for not directly inserting file: to underline later on the file uploaded.
	ifstream input;
	input.open(FileName);

	//Changed data input code so that the name of any text file can be entered by user.
	/*string FileName = "";
	ifstream input;
	cout << "Enter .txt file to open. \n";
	cin >> FileName;
	input.open(FileName);*/

//Initialization of Peptide Data
//----------------------------------------------------------------------------------------------------------------------------------------------//
	cout<<fixed<<setprecision(3);//rounds numbers to second third place.

	int NumberofInput = 0;
	string line;
	/*cout<<"What are the number of inputs ? "<<endl;// Designates array size from number of inputs
	cin>>NumberofInput;*/

	//NumberofInput no longer needs to be manually enterred by user
	//NumberofInputs is found below
	while (getline(input,line))	{
		NumberofInput++;
	}
	cout<<"Number of inputs: " << NumberofInput<<endl;

	input.clear( );
	input.seekg( 0, std::ios::beg );	

	//PatternofBackbones no longer needs to be manually entered by user
	//PatternofBackbones is found below
	int PatternofBackbones = 0;
	string garbage;
	string MolID;
	for (int i = 0; i < NumberofInput; i++)	{
		if (i == 0)	{
			input >> garbage;
			input >> garbage;
			input >> garbage;
			input >> garbage;
			input >> garbage;
			input >> garbage;
		}
		else	{
			PatternofBackbones ++;
			input >> garbage;
			input >> MolID;
			input >> garbage	;
			input >> garbage;
			input >> garbage;
			input >> garbage;
		}
		if (MolID == "BB")	{
			break;
		}
	}

	cout<<"Pattern of backbones: " << PatternofBackbones << endl;

	input.clear( );
	input.seekg( 0, std::ios::beg );
	/*cout<<"Input the pattern of Backbones. Example:Every 'fourth' input is a hydrophilic backbone."<<endl;//Designates array size from (NumberofInputs)/(PatternofBackbones)
	cin>>PatternofBackbones;*/

	string *ResidueID=new string[NumberofInput/PatternofBackbones];//1 or 2; PHE=peptide
	string *MoleculeID=new string[NumberofInput/PatternofBackbones];//Possible Inputs:BB and SC1, SC2, SC3 are discarded
	double *IdentifierNumber=new double[NumberofInput/PatternofBackbones];// input reference number
	double *xcoordinate=new double[NumberofInput/PatternofBackbones];//if you add anything, but a number program shuts down
	double *ycoordinate=new double[NumberofInput/PatternofBackbones];//if you add anything, but a number program shuts down
	double *zcoordinate=new double[NumberofInput/PatternofBackbones];//if you add anything, but a number program shuts down

	//Velocity is irrelevant at this stage
	/*double *xvelocity=new double[NumberofInput/PatternofBackbones];//if you add anything, but a number program shuts down
	double *yvelocity=new double[NumberofInput/PatternofBackbones];//if you add anything, but a number program shuts down
	double *zvelocity=new double[NumberofInput/PatternofBackbones];//if you add anything, but a number program shuts down*/
	int TotalCombinations=(((NumberofInput*NumberofInput)/(PatternofBackbones*PatternofBackbones))/2);//number of combinations for when calculating differences in x, y, and z coordinates
	//TotalCombinations is way more space than needed as seen later on
	double XGreatestValue=0;//greatest length in x
	double YGreatestValue=0;//greatest length in y
	double ZGreatestValue=0;//greatest length in z
	int Placeholder1=0;//placeholder 1 for x,y,z coordinate
	int Placeholder2=0;//placeholder 2 for x,y,z coordinate 
	//double *Xdifference= new double[TotalCombinations];//holds the difference of points in x coordinate 
	//double *Ydifference= new double[TotalCombinations];//holds the difference of points in y coordinate
	//double *Zdifference= new double[TotalCombinations];//holds the difference of points in z coordinate
	double Xaverage=0;//holds average point location in x
	double Yaverage=0;// holds average point location in y
	double Zaverage=0;// holds average point location in z
	int XPlaceholderDifference=0;//placeholder for difference in x coordinate
	int YPlaceholderDifference=0;//placeholder for difference in y coordinate
	int ZPlaceholderDifference=0;//placeholder for difference in z coordinate
	int XOccurences=0;//number of time XGreatestValue occured
	int YOccurences=0;//number of time YGreatestValue occured
	int ZOccurences=0;//number of time zGreatestValue occured
	int CountedInputs=0;// records number of inputs that are hydrophilic inputs
	

//Inputs BB data; non-BB data are sent to the trash which gets rewritten by non-BB data//Number of BBs is counted to match the (NumberofInputs)/(PatternofBackbones) from earlier
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	for(int x=0;x<NumberofInput;x++)	{
		if(x%PatternofBackbones==0)	{
	input>>ResidueID[CountedInputs];
	input>>MoleculeID[CountedInputs];
	input>>IdentifierNumber[CountedInputs];
	input>>xcoordinate[CountedInputs];
	input>>ycoordinate[CountedInputs];
	input>>zcoordinate[CountedInputs];
	/*input>>xvelocity[CountedInputs];
	input>>yvelocity[CountedInputs];
	input>>zvelocity[CountedInputs];*/
	CountedInputs+=1;
	}
	else
	{
	string trash;
	input>>trash;
	input>>trash;
	input>>trash;
	input>>trash;
	input>>trash;
	input>>trash;
	/*input>>trash;
	input>>trash;
	input>>trash;*/
	}
	}

input.close();


//Displays Data that is accounted for in a table format
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	/*cout<<setw(8)<<"Input#"<<setw(11)<<"ResidueId"<<setw(13)<<"Molecule ID"<<setw(14)<<"Identifier #"<<setw(13)<<"xcoordinate"<<setw(13)<<"ycoordinate"<<setw(13)<<"zcoordinate"<<setw(11)<<endl;   
	for(int x=0;x<CountedInputs;x++)
	{   cout<<setw(7)<<left<<x+1<<"|";
	    cout<<setw(10)<<ResidueID[x]<<"|";
		cout<<setw(12)<<left<<MoleculeID[x]<<"|";
		cout<<setw(13)<<left<<fixed<<setprecision(0)<<IdentifierNumber[x]<<fixed<<setprecision(3)<<"|";
		cout<<setw(12)<<left<<xcoordinate[x]<<"|"<<left<<setw(12)<<ycoordinate[x]<<"|"<<setw(12)<<left<<zcoordinate[x]<<setw(10)<<left<<endl;
	}

	cout<<endl;

	cout<<"CountedInputs : "<<CountedInputs<<endl;	*/

//Screens Out Outliers by IQR from arrays, inserts new data set in new arrays and Displays them 
//------------------------------------------------------------------------------------------------------------------------------------------------//
	double *orderedxcoordinate= new double[CountedInputs];
	double *orderedycoordinate= new double[CountedInputs];
	double *orderedzcoordinate= new double[CountedInputs];
	int CurrentMinIndex;
	int n=0;
	int *xplaceholder=new int[CountedInputs];
	int *yplaceholder=new int[CountedInputs];
	int *zplaceholder=new int[CountedInputs];
	double xmedian;
	double xQ1;
	double xQ3;
	double ymedian;
	double yQ1;
	double yQ3;
	double zmedian;
	double zQ1;
	double zQ3;

    //Introduces new array used later to hold original position of values & ordered values
	for(int a=0;a<CountedInputs;a++)
	{
		orderedxcoordinate[a]=xcoordinate[a];
		orderedycoordinate[a]=ycoordinate[a];
		orderedzcoordinate[a]=zcoordinate[a];
	  //cout<<"original xcoordinate contains : "<<originalxcoordinate[a]<<'\t'<<xcoordinate[a]<<'\t'<<"original ycoordinate contains : "<<originalycoordinate[a]<<'\t'<<ycoordinate[a]<<'\t'<<"original zcoordinate contains : "<<originalzcoordinate[a]<<'\t'<<zcoordinate[a]<<'\t'<<a<<endl;
	}
     
	//Puts the xcoordinate into order from greatest to least
	for(int a= 0; a<CountedInputs; a++)
	{
		double CurrentMax = orderedxcoordinate[CountedInputs-1];
		 CurrentMinIndex = CountedInputs-1;

		for(int b=CountedInputs-1; b>=a; b--)
		{ 
			if( CurrentMax< orderedxcoordinate[b])
			{
				CurrentMax = orderedxcoordinate[b];
				CurrentMinIndex = b;
			}
		}

		if (CurrentMinIndex != a)
		{
			orderedxcoordinate[CurrentMinIndex] = orderedxcoordinate[a];	
			orderedxcoordinate[a]= CurrentMax;
		}
	
	}

	//Puts the ycoordinate into order from greatest to least
	for(int a= 0; a<CountedInputs; a++)
	{
		double CurrentMax = orderedycoordinate[CountedInputs-1];
		 CurrentMinIndex = CountedInputs-1;

		for(int b=CountedInputs-1; b>=a; b--)
		{ 
			if( CurrentMax< orderedycoordinate[b])
			{
				CurrentMax = orderedycoordinate[b];
				CurrentMinIndex = b;
			}
		}

		if (CurrentMinIndex != a)
		{
			orderedycoordinate[CurrentMinIndex] = orderedycoordinate[a];	
			orderedycoordinate[a]= CurrentMax;
		}
	
	}

	//Puts the zcoordinate into order from greatest to least
	for(int a= 0; a<CountedInputs; a++)
	{
		double CurrentMax = orderedzcoordinate[CountedInputs-1];
		 CurrentMinIndex = CountedInputs-1;

		for(int b=CountedInputs-1; b>=a; b--)
		{ 
			if( CurrentMax< orderedzcoordinate[b])
			{
				CurrentMax = orderedzcoordinate[b];
				CurrentMinIndex = b;
			}
		}

		if (CurrentMinIndex != a)
		{
			orderedzcoordinate[CurrentMinIndex] = orderedzcoordinate[a];	
			orderedzcoordinate[a]= CurrentMax;
		}
	
	}

	//Displays the new ordered arrays
	/*for (int i = 0; i < CountedInputs; i++)	{
		cout << "Ordered X: " << orderedxcoordinate[i] << endl;
	}
	for (int i = 0; i < CountedInputs; i++)	{
		cout << "Ordered Y: " << orderedycoordinate[i] << endl;
	}
	for (int i = 0; i < CountedInputs; i++)	{
		cout << "Ordered Z: " << orderedzcoordinate[i] << endl;
	}*/

	//Identifies the placeholder of the orderedxcoordinate from original value 
	for(int s=0;s<CountedInputs;s++)
	{
		for(int r=0;r<CountedInputs;r++)
		{
			if(orderedxcoordinate[s]==xcoordinate[r])
			{
				int d=0;
				for(int l=0; l<n;l++)
				{
					if(r==xplaceholder[l])
					{
						d+=1;
					}
				}
				if(d!=0)
				{
					continue;
				}
				xplaceholder[n]=r;
				n+=1;
			}
		}
	}

	//Identifies the placeholder of the orderedycoordinate from original value
	n=0;
	for(int s=0;s<CountedInputs;s++)
	{
		for(int r=0;r<CountedInputs;r++)
		{
			if(orderedycoordinate[s]==ycoordinate[r])
			{
				int d=0;
				for(int l=0; l<n;l++)
				{
					if(r==yplaceholder[l])
					{
						d+=1;
					}
				}
				if(d!=0)
				{
					continue;
				}
				yplaceholder[n]=r;
				n+=1;
			}
		}
	}

	//Identifies the placeholder of the orderedzcoordinate from original value
	n=0;
	for(int s=0;s<CountedInputs;s++)
	{
		for(int r=0;r<CountedInputs;r++)
		{
			if(orderedzcoordinate[r]==zcoordinate[s])
			{
				int d=0;
				for(int l=0; l<n;l++)
				{
					if(r==zplaceholder[l])
					{
						d+=1;
					}
				}
				if(d!=0)
				{
					continue;
				}
				zplaceholder[n]=r;
				n+=1;
			}
		}
	}

	//To check the content
	/*for(int s=0;s<CountedInputs;s++)
	{
		cout<<"ordered x coordinate: "<<orderedxcoordinate[s]<<" placeholder: "<<xplaceholder[s]<<" ordered y coordinate: "<<orderedycoordinate[s]<<" placeholder: "<<yplaceholder[s]<<" ordered z coordinate: "<<orderedzcoordinate[s]<<" placeholder: "<<zplaceholder[s]<<endl;
	}*/

	cout<<'\n'<<endl;
	
//Finds the median of the ordered data set whether even or odd number of data set (seen for the x, y and z)
if((CountedInputs % 2)==0)
{
	int medianplaceholder1=CountedInputs/2;
	int medianplaceholder2=medianplaceholder1+1;
	xmedian=(orderedxcoordinate[medianplaceholder1-1]+orderedxcoordinate[medianplaceholder2-1])/2.0;

	//Two options if odd number of data set and even number of data set for upper and lower area of the data set
	if(((CountedInputs/2) %2 )==0)
	{
	medianplaceholder1=(CountedInputs/2+CountedInputs/4);
	medianplaceholder2=medianplaceholder1+1;
	xQ1=(orderedxcoordinate[medianplaceholder1-1]+orderedxcoordinate[medianplaceholder2-1])/2.0;

	medianplaceholder1=(CountedInputs/2-CountedInputs/4);
	medianplaceholder2=medianplaceholder1+1;
	xQ3=(orderedxcoordinate[medianplaceholder1-1]+orderedxcoordinate[medianplaceholder2-1])/2.0;
	}

	else
	{
		int l=CountedInputs/4;
		int x=((CountedInputs/2))+((CountedInputs/4));
		int y=(CountedInputs/2)-(CountedInputs/4)-1;
		xQ1=orderedxcoordinate[x];
		xQ3=orderedxcoordinate[y];
	}
	
}
else
{
	xmedian=orderedxcoordinate[CountedInputs/2];

	if(((CountedInputs/2) % 2)==0)
	{
		int medianplaceholder1=(CountedInputs/2+1+CountedInputs/4);
		int medianplaceholder2=medianplaceholder1+1;
		xQ1=(orderedxcoordinate[medianplaceholder1-1]+orderedxcoordinate[medianplaceholder2-1])/2.0;

		medianplaceholder1=(CountedInputs/2+1-CountedInputs/4);
		medianplaceholder2=medianplaceholder1-1;
		xQ3=(orderedxcoordinate[medianplaceholder1-1]+orderedxcoordinate[medianplaceholder2-1])/2.0;
	}

	else
	{
		
		xQ1=orderedxcoordinate[((CountedInputs/2)+1)+((CountedInputs/4)+1)-1];
		xQ3=orderedxcoordinate[((CountedInputs/2)+1)-((CountedInputs/4)+1)-1];

	}

}

cout<<"X Median : "<<xmedian<<endl;


if((CountedInputs % 2)==0)
{
	int medianplaceholder1=CountedInputs/2;
	int medianplaceholder2=medianplaceholder1+1;
	ymedian=(orderedycoordinate[medianplaceholder1-1]+orderedycoordinate[medianplaceholder2-1])/2.0;

	if(((CountedInputs/2) %2 )==0)
	{
	medianplaceholder1=(CountedInputs/2+CountedInputs/4);
	medianplaceholder2=medianplaceholder1+1;
	yQ1=(orderedycoordinate[medianplaceholder1-1]+orderedycoordinate[medianplaceholder2-1])/2.0;

	medianplaceholder1=(CountedInputs/2-CountedInputs/4);
	medianplaceholder2=medianplaceholder1+1;
	yQ3=(orderedycoordinate[medianplaceholder1-1]+orderedycoordinate[medianplaceholder2-1])/2.0;
	}

	else
	{
		int l=CountedInputs/4;
		int x=((CountedInputs/2))+((CountedInputs/4));
		int y=(CountedInputs/2)-(CountedInputs/4)-1;
		yQ1=orderedycoordinate[x];
		yQ3=orderedycoordinate[y];
	}
	
}
else
{
	ymedian=orderedycoordinate[CountedInputs/2];

	if(((CountedInputs/2) % 2)==0)
	{
		int medianplaceholder1=(CountedInputs/2+1+CountedInputs/4);
		int medianplaceholder2=medianplaceholder1+1;
		yQ1=(orderedycoordinate[medianplaceholder1-1]+orderedycoordinate[medianplaceholder2-1])/2.0;

		medianplaceholder1=(CountedInputs/2+1-CountedInputs/4);
		medianplaceholder2=medianplaceholder1-1;
		yQ3=(orderedycoordinate[medianplaceholder1-1]+orderedycoordinate[medianplaceholder2-1])/2.0;
	}

	else
	{
		
		yQ1=orderedycoordinate[((CountedInputs/2)+1)+((CountedInputs/4)+1)-1];
		yQ3=orderedycoordinate[((CountedInputs/2)+1)-((CountedInputs/4)+1)-1];

	}

}


cout<<"Y Median : "<<ymedian<<endl;


if((CountedInputs % 2)==0)
{
	int medianplaceholder1=CountedInputs/2;
	int medianplaceholder2=medianplaceholder1+1;
	zmedian=(orderedzcoordinate[medianplaceholder1-1]+orderedzcoordinate[medianplaceholder2-1])/2.0;

	if(((CountedInputs/2) %2 )==0)
	{
	medianplaceholder1=(CountedInputs/2+CountedInputs/4);
	medianplaceholder2=medianplaceholder1+1;
	zQ1=(orderedzcoordinate[medianplaceholder1-1]+orderedzcoordinate[medianplaceholder2-1])/2.0;

	medianplaceholder1=(CountedInputs/2-CountedInputs/4);
	medianplaceholder2=medianplaceholder1+1;
	zQ3=(orderedzcoordinate[medianplaceholder1-1]+orderedzcoordinate[medianplaceholder2-1])/2.0;
	}

	else
	{
		int l=CountedInputs/4;
		int x=((CountedInputs/2))+((CountedInputs/4));
		int y=(CountedInputs/2)-(CountedInputs/4)-1;
		zQ1=orderedzcoordinate[x];
		zQ3=orderedzcoordinate[y];
	}
	
}
else
{
	zmedian=orderedzcoordinate[CountedInputs/2];

	if(((CountedInputs/2) % 2)==0)
	{
		int medianplaceholder1=(CountedInputs/2+1+CountedInputs/4);
		int medianplaceholder2=medianplaceholder1+1;
		zQ1=(orderedzcoordinate[medianplaceholder1-1]+orderedzcoordinate[medianplaceholder2-1])/2.0;

		medianplaceholder1=(CountedInputs/2+1-CountedInputs/4);
		medianplaceholder2=medianplaceholder1-1;
		zQ3=(orderedzcoordinate[medianplaceholder1-1]+orderedzcoordinate[medianplaceholder2-1])/2.0;
	}

	else
	{
		
		zQ1=orderedzcoordinate[((CountedInputs/2)+1)+((CountedInputs/4)+1)-1];
		zQ3=orderedzcoordinate[((CountedInputs/2)+1)-((CountedInputs/4)+1)-1];

	}

}


cout<<"Z Median : "<<zmedian<<endl;
cout<<'\n'<<endl;

double xIQR=xQ3-xQ1;
cout<<"xQ3 : "<<xQ3<<endl;
cout<<"xQ1 : "<<xQ1<<endl;
cout<<"xIQR : "<<xIQR<<endl;

double xmax=xQ3+1.5*xIQR;
double xmin=xQ1-1.5*xIQR;
cout<<"xMax (allowed) : "<<xmax<<endl;
cout<<"xMin (allowed) : "<<xmin<<endl;
int xOutlierCounter=0;
double *xOutliers=new double[CountedInputs];
double *xOutlierPlaceholder=new double[CountedInputs];
double *generalOutlierPlaceholder=new double[CountedInputs];
int generalOutlierCounter=0;
int indicesOfOutliers[CountedInputs];
int curIndexofOutliersIndices = -1;

//initiates the Outlier array and generalOutlierPlaceholder array
for(int s=0;s<CountedInputs;s++)
{
	xOutliers[s]=0.0;//0.0 signifies errors/unused space
	xOutlierPlaceholder[s]=0.0;//0.0 signifies error/unused space
	generalOutlierPlaceholder[s]=-15;//-15 signifies error/unused space
}

//screens out the outliers 
int j=0;//repeat generalOutlierPlaceholderCounter
for(int s=0;s<CountedInputs;s++)
{
	if((xcoordinate[s]>=xmax)||(xcoordinate[s]<=xmin))
	{
		cout<<xcoordinate[s]<<endl;
		cout<<s<<endl;
		curIndexofOutliersIndices++;
		indicesOfOutliers[curIndexofOutliersIndices] = s;
	
		//cout<<"Outlier : "<<orderedxcoordinate[s]<<'\n'<<"Input Number : "<<xplaceholder[s]+1<<endl;
		xOutlierPlaceholder[xOutlierCounter]=xplaceholder[s];
		for(int l=0;l<generalOutlierCounter;l++)
		{
			if(generalOutlierPlaceholder[l]==xplaceholder[s])
			{
				j+=1;
			}
		}
		if(j==0)
		{
			generalOutlierPlaceholder[generalOutlierCounter]=xplaceholder[s];
			generalOutlierCounter+=1;
		}
		xOutliers[xOutlierCounter]=xcoordinate[s];
		xcoordinate[s]=-1;	
		xOutlierCounter+=1;
	}
}

cout << "X Outliers: " << endl;
for (int i = 0; i< xOutlierCounter; i++)	{
	cout << xOutliers[i] << endl;
}

if(xOutlierCounter==0)
{
	cout<<"There are no outliers in xcoordinate."<<endl;
}
else	{
	cout << "There are " << xOutlierCounter <<" outliers in xcoordinate." << endl;
}


int xnonOutlierCounter=0;
for(int s=0; s<CountedInputs;s++)//orderedxcoordinate[s] now doesn't hold any outliers &size = (xnonOutlierCounter) 
{
	if(orderedxcoordinate[s]!=-1)
	{
		//orderedxcoordinate[xnonOutlierCounter]=orderedxcoordinate[s];
		xnonOutlierCounter++;
	}
}

cout<<'\n'<<endl;

double yIQR=yQ3-yQ1;
cout<<"yQ3 : "<<yQ3<<endl;
cout<<"yQ1 : "<<yQ1<<endl;
cout<<"yIQR : "<<yIQR<<endl;

double ymax=yQ3+1.5*yIQR;
double ymin=yQ1-1.5*yIQR;
cout<<"yMax (allowed) : "<<ymax<<endl;
cout<<"yMin (allowed) : "<<ymin<<endl;
int yOutlierCounter=0;
double *yOutliers=new double[CountedInputs];
double *yOutlierPlaceholder=new double[CountedInputs];

for(int s=0;s<CountedInputs;s++)
{
	yOutliers[s]=0.0;//0.0 signifies errors/unused space
	yOutlierPlaceholder[s]=0.0;//0.0 signifies error/unused space
}
j=0;
for(int s=0;s<CountedInputs;s++)
{
	
	if((ycoordinate[s]>=ymax)||(ycoordinate[s]<=ymin))
	{
		bool alreadyThere = false; 
		for (int q = 0; q<curIndexofOutliersIndices+1; q++)	{
			if (indicesOfOutliers[q] == s)	{
				alreadyThere = true;
			}
		}

		if (!alreadyThere)	{
			curIndexofOutliersIndices++;
			indicesOfOutliers[curIndexofOutliersIndices] = s;
			generalOutlierCounter++;
		}

		//cout<<"Outlier : "<<orderedycoordinate[s]<<'\n'<<"Input Number : "<<yplaceholder[s]+1<<endl;
		yOutlierPlaceholder[yOutlierCounter]=yplaceholder[s];
		for(int l=0;l<generalOutlierCounter;l++)
		{
			if(generalOutlierPlaceholder[l]==yplaceholder[s])
			{
				j+=1;
			}
		}
		if(j==0)
		{
			generalOutlierPlaceholder[generalOutlierCounter]=yplaceholder[s];
			//generalOutlierCounter+=1;
		}
		yOutliers[yOutlierCounter]=ycoordinate[s];
		ycoordinate[s]=-1;	
		yOutlierCounter+=1;
	}
}

cout << "Y Outliers: " << endl;
for (int i = 0; i< yOutlierCounter; i++)	{
	cout << yOutliers[i] << endl;
}

if(yOutlierCounter==0)
{
	cout<<"There are no outliers in ycoordinate."<<endl;
}
else	{
	cout << "There are " << yOutlierCounter <<" outliers in ycoordinate." << endl;
}


int ynonOutlierCounter=0;
for(int s=0; s<CountedInputs;s++)
{
	if(ycoordinate[s]!=-1)
	{
		//orderedycoordinate[ynonOutlierCounter]=orderedycoordinate[s];
		ynonOutlierCounter++;
	}
}

cout<<'\n'<<endl;

double zIQR=zQ3-zQ1;
cout<<"zQ3 : "<<zQ3<<endl;
cout<<"zQ1 : "<<zQ1<<endl;
cout<<"zIQR : "<<zIQR<<endl;

double zmax=zQ3+1.5*zIQR;
double zmin=zQ1-1.5*zIQR;
cout<<"zMax (allowed) : "<<zmax<<endl;
cout<<"zMin (allowed) : "<<zmin<<endl;
int zOutlierCounter=0;
double *zOutliers=new double[CountedInputs];
double *zOutlierPlaceholder=new double[CountedInputs];

for(int s=0;s<CountedInputs;s++)
{
	zOutliers[s]=0.0;//0.0 signifies errors/unused space
	zOutlierPlaceholder[s]=0.0;//0.0 signifies error/unused space
}

j=0;
for(int s=0;s<CountedInputs;s++)
{
	
	if((zcoordinate[s]>=zmax)||(zcoordinate[s]<=zmin))
	{

		bool alreadyThere = false; 
		for (int q = 0; q<curIndexofOutliersIndices+1; q++)	{
			if (indicesOfOutliers[q] == s)	{
				alreadyThere = true;
			}
		}

		if (!alreadyThere)	{
			curIndexofOutliersIndices++;
			indicesOfOutliers[curIndexofOutliersIndices] = s;
			generalOutlierCounter++;
		}

		//cout<<"Outlier : "<<orderedzcoordinate[s]<<'\n'<<"Input Number : "<<zplaceholder[s]+1<<endl;
		zOutlierPlaceholder[zOutlierCounter]=zplaceholder[s];
		for(int l=0;l<generalOutlierCounter;l++)
		{
			if(generalOutlierPlaceholder[l]==zplaceholder[s])
			{
				j+=1;
			}
		}
		if(j==0)
		{
			generalOutlierPlaceholder[generalOutlierCounter]=zplaceholder[s];
			//generalOutlierCounter+=1;
		}
		zOutliers[zOutlierCounter]=zcoordinate[s];
		zcoordinate[s]=-1;	
		zOutlierCounter+=1;
	}
}

cout << "Z Outliers: " << endl;
for (int i = 0; i< zOutlierCounter; i++)	{
	cout << zOutliers[i] << endl;
}

if(zOutlierCounter==0)
{
	cout<<"There are no outliers in zcoordinate."<<endl;
}
else	{
	cout << "There are " << zOutlierCounter <<" outliers in zcoordinate." << endl;
}


int znonOutlierCounter=0;
for(int s=0; s<CountedInputs;s++)
{
	if(zcoordinate[s]!=-1)
	{
		//orderedzcoordinate[znonOutlierCounter]=orderedzcoordinate[s];
		znonOutlierCounter++;
	}
}

cout<<'\n'<<endl;

cout<< "Indices of Unique Outliers:"<<endl;
for (int qq = 0; qq < curIndexofOutliersIndices+1; qq++)	{
	//cout<<qq<<endl;
	cout<<indicesOfOutliers[qq]<<endl;
}

//Takes out the outliers from other arrays 
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
/*for(int s=0;s<generalOutlierCounter;s++)
{
	//cout<<"The general outlier placeholders are : "<<generalOutlierPlaceholder[s]<<endl;
}*/

//Assigns generalOutlierCounter to the maximum number of outliers in either x, y, or z
//generalOutlierCounter = std::max(xOutlierCounter,yOutlierCounter);
//generalOutlierCounter = std::max(generalOutlierCounter,zOutlierCounter);
int newCountedInputs=CountedInputs-generalOutlierCounter;
cout<< "Number of new inputs" << newCountedInputs << endl;
cout<<'\n'<<endl;

	string* newResidueID = new string[newCountedInputs];
	string* newMoleculeID = new string[newCountedInputs];
	double* newIdentifierNumber = new double[newCountedInputs];
	double* newxcoordinate = new double[newCountedInputs];
	double* newycoordinate = new double[newCountedInputs];
	double* newzcoordinate = new double[newCountedInputs];
	/*double *newxvelocity=new double[newCountedInputs];
	double *newyvelocity=new double[newCountedInputs];
	double *newzvelocity=new double[newCountedInputs];*/

	std::sort(indicesOfOutliers,indicesOfOutliers+curIndexofOutliersIndices+1);

	j = 0;

	for (int i = 0; i < newCountedInputs;i++)	{
		//cout<<"Current I:"<<endl;
		//cout<<i<<endl;
		while (j<curIndexofOutliersIndices+1)	{
			if (i+j == indicesOfOutliers[j])	{
				j++;
			}
			else	{
				break;
			}
			cout<<endl;
		}

		newResidueID[i] = ResidueID[i+j];
		newMoleculeID[i]=MoleculeID[i+j];
		newIdentifierNumber[i]=IdentifierNumber[i+j];
		newxcoordinate[i]=xcoordinate[i+j];
		newycoordinate[i]=ycoordinate[i+j];
		newzcoordinate[i]=zcoordinate[i+j];

		//j++;
	}

delete[]xcoordinate;
delete[]ycoordinate;
delete[]zcoordinate;
delete[]IdentifierNumber;
delete[]ResidueID;
delete[]MoleculeID;
delete[]orderedxcoordinate;
delete[]orderedycoordinate;
delete[]orderedzcoordinate;
delete[]xplaceholder;
delete[]yplaceholder;
delete[]zplaceholder;
delete[]generalOutlierPlaceholder;


//Displays Outlier-free Data Set
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
/*cout<<"Outlier-Free Data Set (non ordered)"<<endl;
cout<<setw(8)<<"Input#"<<setw(11)<<"ResidueId"<<setw(13)<<"Molecule ID"<<setw(14)<<"Identifier #"<<setw(13)<<"xcoordinate"<<setw(13)<<"ycoordinate"<<setw(13)<<"zcoordinate"<<setw(11)<<endl;   
	for(int x=0;x<newCountedInputs;x++)
	{   cout<<setw(7)<<left<<x+1<<"|";
	    cout<<setw(10)<<newResidueID[x]<<"|";
		cout<<setw(12)<<left<<newMoleculeID[x]<<"|";
		cout<<setw(13)<<left<<fixed<<setprecision(0)<<newIdentifierNumber[x]<<fixed<<setprecision(3)<<"|";
		cout<<setw(12)<<left<<newxcoordinate[x]<<"|"<<left<<setw(12)<<newycoordinate[x]<<"|"<<setw(12)<<left<<newzcoordinate[x]<<setw(10)<<left<<endl;
	}
	
	cout<<'\n'<<"CountedInputs-general Outlier Counter = new Counted Inputs: "<<CountedInputs<<" - "<<generalOutlierCounter<<" = "<<newCountedInputs<<endl;	
	cout<<'\n'<<endl;*/

/*//Computes average location of points=center of structure
//-----------------------------------------------------------------------------------------------------------------------------------------------//
		for(int x=0;x<newCountedInputs;x++)
	{
		Xaverage+=newxcoordinate[x];
		Yaverage+=newycoordinate[x];
		Zaverage+=newzcoordinate[x];
	}
	Xaverage/=newCountedInputs;
	Yaverage/=newCountedInputs;
	Zaverage/=newCountedInputs;

	cout<<"The average x coordinate is : "<<Xaverage<<"."<<endl;
	cout<<"The average y coordinate is : "<<Yaverage<<"."<<endl;
	cout<<"The average z coordinate is : "<<Zaverage<<"."<<endl;
	cout<<'\n'<<endl;*/

//Finds difference of points in (length, width, height)
//----------------------------------------------------------------------------------------------------------------------------------------------//

//This upcoming commented line can be uncommented to see labels of what numbers is getting subtracted 	
	/*cout<<setw(15)<<left<<"x1"<<setw(14)<<"x2"<<setw(14)<<"Difference"<<endl;*/
/*
	while((Placeholder1<(newCountedInputs)))//less than Number of Inputs due to arrays starting at 0 and to not account for Pattern of Backbones causing Placeholder>NumberofInputs
	{   
		if(Placeholder1!=Placeholder2)
		{
			cout<<fixed<<setprecision(0);
			if((-.5<=(newycoordinate[Placeholder1]-newycoordinate[Placeholder2])&&(newycoordinate[Placeholder1]-newycoordinate[Placeholder2])<=.5)&&(-.5<=(newzcoordinate[Placeholder1]-newzcoordinate[Placeholder2])&&(newzcoordinate[Placeholder1]-newzcoordinate[Placeholder2])<=.5))
			//The above .2 and .5 can be changed as this is the amount of acceptable variation in y & z to be able to get a length in x
			{
				cout<<fixed<<setprecision(3);
				Xdifference[XPlaceholderDifference]=newxcoordinate[Placeholder1]-newxcoordinate[Placeholder2];// Finds the difference in x coordinates

				if(Xdifference[XPlaceholderDifference]>0)
				{
				cout<<setw(15)<<left<<newxcoordinate[Placeholder1]<<setw(14)<<left<<newxcoordinate[Placeholder2]<<setw(14)<<Xdifference[XPlaceholderDifference]<<endl;
				XPlaceholderDifference=XPlaceholderDifference+1;
				}
			}
		}
		if(Placeholder2==(newCountedInputs-1))
	    {
			Placeholder1=Placeholder1+1;
	        Placeholder2=0;
	    }
	     
		else
	    {
		Placeholder2=Placeholder2+1;
	    }
	      
	}

	cout<<'\n'<<endl;

	Placeholder1=0;
	Placeholder2=0;

	cout<<setw(15)<<left<<"y1"<<setw(14)<<"y2"<<setw(14)<<"Difference"<<endl;
	while((Placeholder1<(newCountedInputs)))//less than Number of Inputs due to arrays starting at 0 and to not account for Pattern of Backbones causing Placeholder>NumberofInputs
	{   
		if(Placeholder1!=Placeholder2)
		{
			cout<<fixed<<setprecision(0);
			if((-.5<=(newxcoordinate[Placeholder1]-newxcoordinate[Placeholder2])&&(newxcoordinate[Placeholder1]-newxcoordinate[Placeholder2])<=.5)&&(-.5<=(newzcoordinate[Placeholder1]-newzcoordinate[Placeholder2])&&(newzcoordinate[Placeholder1]-newzcoordinate[Placeholder2])<=.5))
			{
				cout<<fixed<<setprecision(3);
				Ydifference[YPlaceholderDifference]=newycoordinate[Placeholder1]-newycoordinate[Placeholder2];// Finds the difference in y coordinates
			    
				if(Ydifference[YPlaceholderDifference]>0)
				{
				cout<<setw(15)<<left<<newycoordinate[Placeholder1]<<setw(14)<<left<<newycoordinate[Placeholder2]<<setw(14)<<Ydifference[YPlaceholderDifference]<<endl;
				YPlaceholderDifference=YPlaceholderDifference+1;
				}
			}
		}
		if(Placeholder2==(newCountedInputs-1))
	    {
			Placeholder1=Placeholder1+1;
	        Placeholder2=0;
	    }
	     
		else
	    {
		Placeholder2=Placeholder2+1;
	    }
	      
	}

	cout<<'\n'<<endl;

	Placeholder1=0;
	Placeholder2=0;
	
	cout<<setw(15)<<left<<"z1"<<setw(14)<<"z2"<<setw(14)<<"Difference"<<endl;
	while((Placeholder1<(newCountedInputs)))//less than Number of Inputs due to arrays starting at 0 and to not account for Pattern of Backbones causing Placeholder>NumberofInputs
	{   
		if(Placeholder1!=Placeholder2)
		{
			cout<<fixed<<setprecision(0);
			if((-.5<=(newxcoordinate[Placeholder1]-newxcoordinate[Placeholder2])&&(newxcoordinate[Placeholder1]-newxcoordinate[Placeholder2])<=.5)&&(-.5<=(newycoordinate[Placeholder1]-newycoordinate[Placeholder2])&&(newycoordinate[Placeholder1]-newycoordinate[Placeholder2])<=.5))
			{
			    cout<<fixed<<setprecision(3);
				Zdifference[ZPlaceholderDifference]=newzcoordinate[Placeholder1]-newzcoordinate[Placeholder2];// Finds the difference in z coordinates
				
				if(Zdifference[ZPlaceholderDifference]>0)
				{
				cout<<setw(15)<<left<<newzcoordinate[Placeholder1]<<setw(14)<<left<<newzcoordinate[Placeholder2]<<setw(14)<<Zdifference[ZPlaceholderDifference]<<endl;
				ZPlaceholderDifference=ZPlaceholderDifference+1;
				}
			}
		}
		if(Placeholder2==(newCountedInputs-1))
	    {
			Placeholder1=Placeholder1+1;
	        Placeholder2=0;
	    }
	     
		else
	    {
		Placeholder2=Placeholder2+1;
	    }
	      
	}

	cout<<'\n'<<endl;
	cout<<fixed<<setprecision(3);
	*/

//Finds the average length in x,y, and z
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	/*double sumlengthx=0;
	double AverageLengthx;
	for(int r=0;r<XPlaceholderDifference;r++)
	{
		sumlengthx+=Xdifference[r];
	}
	AverageLengthx=sumlengthx/XPlaceholderDifference;

	double sumlengthy=0;
	double AverageLengthy;
	for(int r=0;r<YPlaceholderDifference;r++)
	{
		
		sumlengthy+=Ydifference[r];
	}
	AverageLengthy=sumlengthy/YPlaceholderDifference;

	double sumlengthz=0;
	double AverageLengthz;
	for(int r=0;r<ZPlaceholderDifference;r++)
	{
		sumlengthz+=Zdifference[r];
	}
	AverageLengthz=sumlengthz/ZPlaceholderDifference;

	cout<<"The average length in x : "<<AverageLengthx<<" in y : "<<AverageLengthy<<" in z : "<<AverageLengthz<<"."<<endl;

//Finds the greatest length in x, y, and z
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	double XGreatestLength=0;
	int XGreatestLengthPlaceholder; // why does this matter?
	for(int r=0;r<XPlaceholderDifference;r++)
	{
		if(XGreatestLength<Xdifference[r])
		{
			XGreatestLength=Xdifference[r];
			XGreatestLengthPlaceholder=r;
		}
	}

	double YGreatestLength=0;
	int YGreatestLengthPlaceholder; //why does this matter?
	for(int r=0;r<YPlaceholderDifference;r++)
	{
		if(YGreatestLength<Ydifference[r])
		{
			YGreatestLength=Ydifference[r];
			YGreatestLengthPlaceholder=r;
		}
	}

	double ZGreatestLength=0;
	int ZGreatestLengthPlaceholder; //why does this matter?
	for(int r=0;r<ZPlaceholderDifference;r++)
	{
		if(ZGreatestLength<Zdifference[r])
		{
			ZGreatestLength=Zdifference[r];
			ZGreatestLengthPlaceholder=r;
		}
	}
	
	cout<<"The greatest length in x : "<<XGreatestLength<<" in y : "<<YGreatestLength<<" in z : "<<ZGreatestLength<<"."<<endl;

//Find the least length in x,y, and z
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
	double XLeastLength=XGreatestLength;
	int XLeastLengthPlaceholder; // why does this matter?
	for(int r=0;r<XPlaceholderDifference;r++)
	{
		if(XLeastLength>Xdifference[r])
		{
			XLeastLength=Xdifference[r];
			XLeastLengthPlaceholder=r;
		}
	}

	double YLeastLength=YGreatestLength;
	int YLeastLengthPlaceholder; //why does this matter?
	for(int r=0;r<YPlaceholderDifference;r++)
	{
		if(YLeastLength>Ydifference[r])
		{
			YLeastLength=Ydifference[r];
			YLeastLengthPlaceholder=r;
		}
	}

	double ZLeastLength=ZGreatestLength;
	int ZLeastLengthPlaceholder; //why does this matter?
	for(int r=0;r<ZPlaceholderDifference;r++)
	{
		if(ZLeastLength>Zdifference[r])
		{
			ZLeastLength=Zdifference[r];
			ZLeastLengthPlaceholder=r;
		}
	}
	
	cout<<"The least length in x : "<<XLeastLength<<" in y : "<<YLeastLength<<" in z : "<<ZLeastLength<<"."<<endl;

	cout<<'\n'<<endl;
	
// orient the structure's points to a center at 0,0,0
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
for(int r=0;r<newCountedInputs;r++)
{
	newxcoordinate[r]=newxcoordinate[r]-Xaverage;
	newycoordinate[r]=newycoordinate[r]-Yaverage;
	newzcoordinate[r]=newzcoordinate[r]-Zaverage;
}*/

//Displays oriented Data
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
/*cout<<setw(8)<<"Input#"<<setw(11)<<"ResidueId"<<setw(13)<<"Molecule ID"<<setw(14)<<"Identifier #"<<setw(13)<<"xcoordinate"<<setw(13)<<"ycoordinate"<<setw(13)<<"zcoordinate"<<setw(11)<<"xvelocity"<<setw(11)<<"yvelocity"<<setw(11)<<"zvelocity"<<endl;   
	for(int x=0;x<newCountedInputs;x++)
	{   cout<<setw(7)<<left<<x+1<<"|";
	    cout<<setw(10)<<newResidueID[x]<<"|";
		cout<<setw(12)<<left<<newMoleculeID[x]<<"|";
		cout<<setw(13)<<left<<fixed<<setprecision(0)<<newIdentifierNumber[x]<<fixed<<setprecision(3)<<"|";
		cout<<setw(12)<<left<<newxcoordinate[x]<<"|"<<left<<setw(12)<<newycoordinate[x]<<"|"<<setw(12)<<left<<newzcoordinate[x]<<"|"<<setw(10)<<left<<newxvelocity[x]<<"|"<<left<<setw(10)<<newyvelocity[x]<<"|"<<setw(10)<<left<<newzvelocity[x]<<"|"<<endl;
	}
	cout<<endl;*/

//Find center of mass (assuming mass of each particle is 1)
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
/*
double xcm;
double ycm;
double zcm;
double xsum=0;
double ysum=0;
double zsum=0;

for(int r=0; r<newCountedInputs;r++)
{
	xsum+=newxcoordinate[r];
	ysum+=newycoordinate[r];
	zsum+=newzcoordinate[r];
}
xcm=xsum/static_cast<double>(newCountedInputs);
ycm=ysum/static_cast<double>(newCountedInputs);
zcm=zsum/static_cast<double>(newCountedInputs);

cout<<"The new center of mass is at the origin. x : "<<xcm<<" y : "<<ycm<<" z : "<<zcm<<'\n'<<endl;
*/

//Determining if Nanovesicle
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

//Finding number of water data points
string WaterFileName="solventTextFile.txt"; //Reading .txt file with just water data
ifstream input2;
input2.open(WaterFileName);

int NumberofWaterInput = 0;
string line2;

while (getline(input2,line2))	{
	NumberofWaterInput++;
}

//Subtract one to account for box dimensions line at the end
NumberofWaterInput = NumberofWaterInput - 1;
cout<<"Number of water inputs: " << NumberofWaterInput<<endl;

input2.clear( );
input2.seekg( 0, std::ios::beg );

string trash;

//This first pass through the data is to obtain the box dimensions at bottom of solvent file
for(int x=newCountedInputs; x < newCountedInputs + NumberofWaterInput; x++)	{
	input2>>trash;
	input2>>trash;
	if ((x - newCountedInputs) < (10000 - NumberofInput - 1))	{
		input2>>trash;
	}

	input2>>trash;
	input2>>trash;
	input2>>trash;
}

double* dimensions = new double[3];
for (int i = 0; i < 3; i++)	{
	input2>>dimensions[i];
}

input2.clear( );
input2.seekg( 0, std::ios::beg );

cout << dimensions[0] << endl;
cout << dimensions[1] << endl;
cout << dimensions[2] << endl;

//Move center of system to center of box
double xsum = 0;
double ysum = 0;
double zsum = 0;
double xAverage = 0;
double yAverage = 0;
double zAverage = 0;

for (int i = 0; i < (newCountedInputs); i++)	{
	xsum = xsum + newxcoordinate[i];
	ysum = ysum + newycoordinate[i];
	zsum = zsum + newzcoordinate[i];
}

xAverage = xsum/(newCountedInputs);
yAverage = ysum/(newCountedInputs);
zAverage = zsum/(newCountedInputs);

cout << "x" << xAverage<< endl;
cout << "y" << yAverage<< endl;
cout << "z" << zAverage<< endl;

for(int i = 0; i < (newCountedInputs); i++)	{
	double temp = newxcoordinate[i];
	double temp2 = newycoordinate[i];
	double temp3 = newzcoordinate[i];
	newxcoordinate[i] = (temp-xAverage + dimensions[0]/2);
	newycoordinate[i] = (temp2-yAverage + dimensions[1]/2);
	newzcoordinate[i] = (temp3-zAverage + dimensions[2]/2);
}

//Confirms that new center of sysytem is at center of box
double xsum2 = 0;
double xAverage2;
double ysum2 = 0;
double yAverage2;
double zsum2 = 0;
double zAverage2;

for (int i = 0; i < (newCountedInputs); i++)	{
	xsum2 += newxcoordinate[i];
	ysum2 += newycoordinate[i];
	zsum2 += newzcoordinate[i];
}

xAverage2 = xsum2/(newCountedInputs);
yAverage2 = ysum2/ (newCountedInputs);
zAverage2 = zsum2/ (newCountedInputs);


cout<<"The new center of mass is at the center of box. x : "<<xAverage2<<" y : "<<yAverage2<<" z : "<<zAverage2<<'\n'<<endl;

//Initializing and copying peptide-only data into arrays that will contain peptide and water data
string* ResidueID2 = new string[NumberofWaterInput + newCountedInputs];
string* MoleculeID2 = new string[NumberofWaterInput + newCountedInputs];
double* IdentifierNumber2 = new double[NumberofWaterInput + newCountedInputs];
double* xcoordinate2 = new double[NumberofWaterInput + newCountedInputs];
double* ycoordinate2 = new double[NumberofWaterInput + newCountedInputs];
double* zcoordinate2 = new double[NumberofWaterInput + newCountedInputs];

for (int i = 0; i < newCountedInputs; i++)	{
	ResidueID2[i] = newResidueID[i];
	MoleculeID2[i] = newMoleculeID[i];
	IdentifierNumber2[i] = newIdentifierNumber[i];
	xcoordinate2[i] = newxcoordinate[i];
	ycoordinate2[i] = newycoordinate[i];
	zcoordinate2[i] = newzcoordinate[i];
}

delete[]newResidueID;
delete[]newMoleculeID;
delete[]newIdentifierNumber;
delete[]newxcoordinate;
delete[]newycoordinate;
delete[]newzcoordinate;

/*std::copy(newResidueID, newResidueID + newCountedInputs, ResidueID2);
std::copy(newMoleculeID, newMoleculeID + newCountedInputs,MoleculeID2);
std::copy(newIdentifierNumber, newIdentifierNumber + newCountedInputs, IdentifierNumber2);
std::copy(newxcoordinate, newxcoordinate + newCountedInputs, xcoordinate2);
std::copy(newycoordinate, newycoordinate + newCountedInputs, ycoordinate2);
std::copy(newzcoordinate, newzcoordinate + newCountedInputs, zcoordinate2);*/

//string trash;

//Read in water data
for(int x=newCountedInputs; x < newCountedInputs + NumberofWaterInput; x++)	{
	//input2>>ResidueID2[x];
	//input2>>MoleculeID2[x];
	//input2>>IdentifierNumber2[x];
	input2>>trash;
	input2>>trash;
	if ((x - newCountedInputs) < (10000 - NumberofInput - 1))	{
		input2>>trash;
	}

	input2>>xcoordinate2[x];
	input2>>ycoordinate2[x];
	input2>>zcoordinate2[x];

	ResidueID2[x] = "solvent";
	MoleculeID2[x] = "S";
	IdentifierNumber2[x] = x+1;
}

//Read box dimensions from end of solvent file
/*double* dimensions = new double[3];
for (int i = 0; i < 3; i++)	{
	input2>>dimensions[i];
}*/


input2.close();

//Displays Outlier-free Data Set
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
cout<<endl;
cout<<"Data Set with Water Particles"<<endl;
cout<<setw(8)<<"Input#"<<setw(11)<<"ResidueId"<<setw(13)<<"Molecule ID"<<setw(14)<<"Identifier #"<<setw(13)<<"xcoordinate"<<setw(13)<<"ycoordinate"<<setw(13)<<"zcoordinate"<<setw(11)<<endl;   
/*for(int x=0;x<(newCountedInputs + NumberofWaterInput);x++)
{   cout<<setw(7)<<left<<x+1<<"|";
    cout<<setw(10)<<ResidueID2[x]<<"|";
	cout<<setw(12)<<left<<MoleculeID2[x]<<"|";
	cout<<setw(13)<<left<<fixed<<setprecision(0)<<IdentifierNumber2[x]<<fixed<<setprecision(3)<<"|";
	cout<<setw(12)<<left<<xcoordinate2[x]<<"|"<<left<<setw(12)<<ycoordinate2[x]<<"|"<<setw(12)<<left<<zcoordinate2[x]<<setw(10)<<left<<endl;
}
cout<<"\n"<<endl;*/

//Calls the nanovesicle structure function to check if it's a nanovesicle
//double tempDimensions[3] = {2.0,3.0,4.0}; //Temporary placeholder for the dimensions arrays for now until actual box dimensions can be read from the txt file

if (waterTrappingCheck(xcoordinate2, ycoordinate2, zcoordinate2, ResidueID2, MoleculeID2, IdentifierNumber2, NumberofWaterInput + newCountedInputs, dimensions))	{
	cout << "Structure traps water." << endl;
	ofstream myfile;
	myfile.open ("output.txt");
	if (Nanovesicle(newxcoordinate, newycoordinate, newzcoordinate, newResidueID, newMoleculeID, newIdentifierNumber, newCountedInputs, dimensions))	{
		cout << "Nanovesicle" << endl;
		myfile << "v";
		
	}
	else	{
		cout << "Disordered" << endl;
		myfile << "d";
	}
	myfile.close();
}

//finding (a,0,0)/(0,b,0)/(0,0,c)

//Determining Bilayer Structure
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//Can improve Bilayer by changing the the 3 points used to making the equations as only the first three are used
/*double Point1[3];
double Point2[3];
double Point3[3];
double distance1[3];
double distance2[3];
double normalvector[3];
double BilayerSuccessCounter=0;
double GreatestPoint1[3];
double GreatestPoint2[3];
double GreatestPoint3[3];
double GreatestBilayerSuccessCounter=0;


for(int Point3Placeholder=0;Point3Placeholder<newCountedInputs;Point3Placeholder++)
{
	for(int Point2Placeholder=0;Point2Placeholder<newCountedInputs;Point2Placeholder++)
	{
		for(int Point1Placeholder=0;Point1Placeholder<newCountedInputs;Point1Placeholder++)
		{
			

					Point1[0]=newxcoordinate[Point1Placeholder];
					Point1[1]=newycoordinate[Point1Placeholder];
					Point1[2]=newzcoordinate[Point1Placeholder];
				
					Point2[0]=newxcoordinate[Point2Placeholder];
					Point2[1]=newycoordinate[Point2Placeholder];
					Point2[2]=newzcoordinate[Point2Placeholder];
				
					Point3[0]=newxcoordinate[Point3Placeholder];
					Point3[1]=newycoordinate[Point3Placeholder];
					Point3[2]=newzcoordinate[Point3Placeholder];
				
if(((Point1[0]!=Point2[0])||(Point1[1]!=Point2[1])||(Point1[2]!=Point2[2]))&&((Point1[0]!=Point3[0])||(Point1[1]!=Point3[1])||(Point1[2]!=Point3[2]))&&((Point3[0]!=Point2[0])||(Point3[1]!=Point2[1])||(Point3[2]!=Point2[2])))
{
		for(int r=0;r<3;r++)//Point1 considered the vertex
			{
				distance1[r]=Point2[r]-Point1[r];
				distance2[r]=Point3[r]-Point1[r];
			}

				normalvector[0]=(distance1[1]*distance2[2]-distance1[2]*distance2[1]);
				normalvector[1]=(distance1[0]*distance2[2]-distance1[2]*distance2[0])*-1;
				normalvector[2]=(distance1[0]*distance2[1]-distance1[1]*distance2[0]);


			double equationofbilayer;

			for(int r=0;r<newCountedInputs;r++)//Point1 is used as initial point
			{
				equationofbilayer=normalvector[0]*(newxcoordinate[r]-Point1[0])+normalvector[1]*(newycoordinate[r]-Point1[1])+normalvector[2]*(newzcoordinate[r]-Point1[2]);
				if((-.1<equationofbilayer)&(equationofbilayer<.1))//up to you for acceptable range
					{
						BilayerSuccessCounter+=1;
					}
			}
			if(GreatestBilayerSuccessCounter<BilayerSuccessCounter)
			{
				GreatestBilayerSuccessCounter=BilayerSuccessCounter;
				GreatestPoint1[0]=Point1[0];
				GreatestPoint1[1]=Point1[1];
				GreatestPoint1[2]=Point1[2];
				GreatestPoint2[0]=Point2[0];
				GreatestPoint2[1]=Point2[1];
				GreatestPoint2[2]=Point2[2];
				GreatestPoint3[0]=Point3[0];
				GreatestPoint3[1]=Point3[1];
				GreatestPoint3[2]=Point3[2];
			}
			//cout<<BilayerSuccessCounter/newCountedInputs*100<<"%"<<endl;
			BilayerSuccessCounter=0;
			}
		}
	}	
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
cout<<endl<<'\n';
cout<<"Percentage of Points that fits Nanobilayer structure : "<<GreatestBilayerSuccessCounter/newCountedInputs*100<<" %"<<endl;
cout<<"Greatest Point 1 : "<<GreatestPoint1[0]<<","<<GreatestPoint1[1]<<","<<GreatestPoint1[2]<<endl<<"Greatest Point 2 : "<<GreatestPoint2[0]<<","<<GreatestPoint2[1]<<","<<GreatestPoint2[2]<<endl
	<<"Greatest Point 3 : "<<GreatestPoint3[0]<<","<<GreatestPoint3[1]<<","<<GreatestPoint3[2]<<endl;
if(FileName=="FFBilayer.txt")
{
	cout<<right<<setw(63)<<"________"<<endl<<left;
}

cout<<endl;


//Nanotube
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//Finding the lowest positive number
double a2;
double b2;
double c2;
double a3;
double b3;
double c3;
double lowesty;
double lowestx;
double lowestz;
int lowestyplaceholder;
int lowestxplaceholder;
int lowestzplaceholder;

for(int r=0; r<newCountedInputs;r++)
{
	if(r==0)
	{
		lowestz=newzcoordinate[r];
		lowestzplaceholder=r;
		lowesty=newycoordinate[r];
		lowestyplaceholder=r;
		lowestx=newxcoordinate[r];
		lowestxplaceholder=r;
		
	}
	if((0<=newzcoordinate[r])&((newzcoordinate[r]<lowestz)||(lowestz<0)))
	{
		lowestz=newzcoordinate[r];
		lowestzplaceholder=r;
	}
	if((0<=newycoordinate[r])&((newycoordinate[r]<lowesty)||(lowesty<0)))
	{
		lowesty=newycoordinate[r];
		lowestyplaceholder=r;
	}
	if((0<=newxcoordinate[r])&((newxcoordinate[r]<lowestx)||(lowestx<0)))
	{
		lowestx=newxcoordinate[r];
		lowestxplaceholder=r;
	}
}
a2=newxcoordinate[lowestyplaceholder];
b2=newycoordinate[lowestxplaceholder];

//opening centered on z-axis

double equationofznanotube;
double ZNanotubeSuccessCounter=0;

for(int r=0;r<newCountedInputs;r++)
{
	equationofznanotube=(pow(newxcoordinate[r],2)/pow(a2,2))+(pow(newycoordinate[r],2)/pow(b2,2));
	
	if((.5<equationofznanotube)&(equationofznanotube<1.5))
	{
		ZNanotubeSuccessCounter+=1;
	}
}

cout<<"Percentage of Points that fits Nanotube structure opening on z-axis: "<<ZNanotubeSuccessCounter/newCountedInputs*100<<" %"<<endl;
if(FileName=="FFNanotube.txt")
{
	cout<<right<<setw(77)<<"________"<<endl<<left;
}
cout<<endl;

//opening centered on y-axis


a3=newxcoordinate[lowestzplaceholder];
c2=newzcoordinate[lowestxplaceholder];

double equationofynanotube;
double YNanotubeSuccessCounter=0;

for(int r=0;r<newCountedInputs;r++)
{
	equationofynanotube=(pow(newxcoordinate[r],2)/pow(a3,2))+(pow(newzcoordinate[r],2)/pow(c2,2));
	
	if((.5<equationofynanotube)&(equationofynanotube<1.5))
	{
		YNanotubeSuccessCounter+=1;
	}
}

cout<<"Percentage of Points that fits Nanotube structure opening on y-axis: "<<YNanotubeSuccessCounter/newCountedInputs*100<<" %"<<endl;
if(FileName=="FFNanotube.txt")
{
	cout<<right<<setw(77)<<"________"<<endl<<left;
}

	
cout<<endl;


//opening centered on x-axis


c3=newzcoordinate[lowestyplaceholder];
b3=newycoordinate[lowestzplaceholder];

double equationofxnanotube;
double XNanotubeSuccessCounter=0;

for(int r=0;r<newCountedInputs;r++)
{
	equationofxnanotube=(pow(newzcoordinate[r],2)/pow(c3,2))+(pow(newycoordinate[r],2)/pow(b3,2));
	
	if((.5<equationofxnanotube)&(equationofxnanotube<1.5))
	{
		XNanotubeSuccessCounter+=1;
	}
}

cout<<"Percentage of Points that fits Nanotube structure opening on x-axis: "<<XNanotubeSuccessCounter/newCountedInputs*100<<" %"<<endl;
if(FileName=="FFNanotube.txt")
{
	cout<<right<<setw(77)<<"________"<<endl<<left;
}


	
cout<<endl;*/

//Defining Characterisitics of structures
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//

/*if(XNanotubeSuccessCounter>
*/

//Distribution Table/Graph
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//used for value ranges
/*
int greatestdistancex=0;
int greatestdistancey=0;
int greatestdistancez=0;
//First set the range by finding the length for the distribution graph
for(int r=0;r<newCountedInputs;r++)
{
		if(greatestdistancex<newxcoordinate[r])
		{
			greatestdistancex=newxcoordinate[r];
		}
}
for(int r=0;r<newCountedInputs;r++)
{
		if(greatestdistancey<newycoordinate[r])
		{
			greatestdistancey=newycoordinate[r];
		}
}
for(int r=0;r<newCountedInputs;r++)
{
		if(greatestdistancez<newzcoordinate[r])
		{
			greatestdistancez=newzcoordinate[r];
		}
}
int xdistributionarraycounter=0;
for(int r=-greatestdistancex;r<=greatestdistancex;r++)
{
	xdistributionarraycounter++;
}
//actual data collection
int *xdistributionarray=new int[xdistributionarraycounter];
for(int r=-greatestdistancex;r<=greatestdistancex;r++)
{
	for(int s=0;s<newCountedInputs;s++)//x view collection
	{
		if((newxcoordinate[s]>=(-r+.02)) && (newxcoordinate[s]<=(-r-.02)))
		{
			if((newycoordinate[s]>=-.02) &&(newycoordinate[s]<=.02))
			{
			}
		}
}
	*/

	//input.close();
	/*delete[]ResidueID2;
	delete[]MoleculeID2;
	delete[]IdentifierNumber2;
	delete[]MoleculeID2;
	delete[]IdentifierNumber2;
	delete[]xcoordinate2;
	delete[]ycoordinate2;
	delete[]zcoordinate2;*/

	return 0;
}