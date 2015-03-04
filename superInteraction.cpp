#include <iostream>        //för cout m.m
#include <cstdlib>         //för abort()-funktionen.
#include <vector>          //för vektorklassen
#include <cmath>           //för pow-, sqrt-funktionerna
#include <cassert>         //for all asserts()

#include "auxiliary.h"
#include "superInteraction.h"

using namespace std;

//for testing
const bool VERBOSE = 0;

//NOTE: This code will not work for fixed boundary (yet).

//constructor for the SuperInteraction-class initiate the constructor
//of the lattice-class:
SuperInteraction::SuperInteraction(int xSquare,int ySquare,int zSquare,
                                   int particleNumber, double seed,
                                   bool boundary,float interaction)
: Lattice(xSquare,ySquare,zSquare, particleNumber, seed, boundary){

  //other specific stuff for my SuperInteraction-code.
  interactionStrength_ = interaction;

}

void SuperInteraction::setInteraction(float interaction) const{
  //override the lattice-version, since using that version
  //would be fatal. This one doesn't do anything, other than
  //overriding the former.
  cerr << "Not using setInteraction()-function"
       << "interaction strength is set in constructor" << endl;
}



void SuperInteraction::printPositionToScreen(int particle){
  //Print position of particle label (int), or vector
  //containing particle labels, to screen. (Test-function)
  int m;
  m = particle;
  cout << "Particle " << m << " X:Y:Z " << pos_[m] << endl;
}

void SuperInteraction::printPositionToScreen(vector<int> particle){
  //Print position of particle label (int), or vector
  //containing particle labels, to screen. (Test-function)
  int m;
  for(size_t i = 0; i < particle.size(); i++){
    m = particle[i];
    cout << "Particle " << m << " X:Y:Z " << pos_[m] << endl;
  }
  cout << "Total size: " << particle.size() << endl;
}



void SuperInteraction::controlVacancyCheck3(int n){
  //This is just a control/bug finding algorithm, based on VacancyCheck2()
  //It's run before setting the new site as occupied, to make sure it was
  //vacant before the move...

  if(n < noParticles_ && 0 <= n){
    if(board_.isOccupied(pos_[n])){  //if occupied
      cout << endl << "Occupied (" << pos_[n] << ")" << endl;
      throw std::string("Two particles on same site.");
    }
  }
  else
    throw std::string("Accessing invalid particle.");
}


void  SuperInteraction::generateTrajectory(int& E){
  //This is just a wrapper for easier passing to the thread-library.

  place();  //place particles on lattice
  move();   //perform a simulation

  //increase the iterator that counts ensembles:
  ++E;
}


double SuperInteraction::sumFrictionCoefficient(vector<int> cluster) const{
  //Function to sum the friction coefficients (inverse jump rate) of
  //all particles in the cluster.

  double sum = 0;
  int vectorLength = (int) cluster.size();

  for(int i = 0; i < vectorLength; i++){
      int m = cluster[i]; //gives correct particle-labels
      sum += 1.0/jumpRate_[m].x.r;
      sum += 1.0/jumpRate_[m].x.l;
  }

  if(dim_ >= 2){
    for(int i = 0; i < vectorLength; i++){
      int m = cluster[i]; //gives correct particle-labels
      sum += 1.0/jumpRate_[m].y.r;
      sum += 1.0/jumpRate_[m].y.l;
    }
  }

  if(dim_ == 3){
    for(int i = 0; i < vectorLength; i++){
      int m = cluster[i]; //gives correct particle-labels
      sum += 1.0/jumpRate_[m].z.l;
      sum += 1.0/jumpRate_[m].z.r;
    }
  }

  return sum;
}


double SuperInteraction::sumFrictionCoefficient(int particle) const{
  //Function to sum the friction coefficients (inverse jump rate) of
  //all particles in the cluster.

  double sum = 0;
  int vectorLength = 1;
  int m = particle;

  //(Don't actually need the for-loop when just one single particle)
  for(int i = 0; i < vectorLength; i++){
    sum += 1.0/jumpRate_[m].x.r;
    sum += 1.0/jumpRate_[m].x.l;
  }

  if(dim_ >= 2){
    for(int i = 0; i < vectorLength; i++){
      sum += 1.0/jumpRate_[m].y.r;
      sum += 1.0/jumpRate_[m].y.l;
    }
  }

  if(dim_ == 3){
    for(int i = 0; i < vectorLength; i++){
      sum += 1.0/jumpRate_[m].z.l;
      sum += 1.0/jumpRate_[m].z.r;
    }
  }

  return sum;
}






void SuperInteraction::superInteractionCode(int n, int direction){
  //Note: this function should not be called from MoveAndBoundaryCheck, but
  //rather from the function that determines n, namely Move2()

  //This is a switch to print the events to the screen,
  //to check it's actually doing what we want it to do.
  bool verbose = VERBOSE;

  if(verbose){
    string strDir[] = {"X+1","X-1","Y+1","Y-1","Z+1","Z-1"};
    cout << endl << "-------------------------" << endl
         << "New run of SuperInteractionCode() to move particle n = "
         << n <<", in direction: " << strDir[direction] << endl;
    printPositionToScreen(n);
    cout << "-------------------------" << endl << endl;

    //print all occupied sites to screen:
    cout<<"Lattice configuration is:"<<endl;
    for(int i = 0; i < noParticles_; i++){
      cout << "particle " << i << ": " << pos_[i].x <<", " << pos_[i].y;
      if(dim_ == 3) cout <<", "<< pos_[i].z; cout << endl;
    }
    cout << endl;
  }

  assert(board_.checkOK());

  vector<int> cluster;       //store number of particles in Cluster here.
  board_.checkNeighbours(pos_[n], cluster, isBoundaryFix_);
  double V = interactionStrength_;  //interaction energy /(k_B*T) (V > 0)

  if(cluster.empty()){
    //accept move
    board_.setAsVacant(pos_[n]);
    moveAndBoundaryCheck(n,direction);           //move
    controlVacancyCheck3(n);                     //TEST that vacancy(n)..=-1
    board_.setAsOccupied(pos_[n], n);            //new position occupied
    assert(board_.checkOK());

    if(verbose){
      cout << "Moved single particle (no cluster/neighbors) n to " << endl;
      printPositionToScreen(n);
      cout << endl;
    }

  }
  else{
    //start with adding particle n to the cluster:
    //(or else it might not be part of the cluster,
    // due to Boltzmann prob, which would be wrong!)

    //    Cluster.push_back(n);
    //    Gör detta i BuildCluster2() istället

    assert(board_.checkOK());
    buildCluster2(n,cluster,V);

    //TEST---------
    bool found_n = false;
    for(size_t i = 0; i < cluster.size(); i++){
      if(cluster[i] == n)
        found_n = true;
    }
    assert(found_n);
    //-------------


    if(verbose){
      cout << "Particle "<< n <<" part of cluster:" << endl;
      printPositionToScreen(cluster);
      cout << endl;
    }

    //Calculate the friction constant (k^(-1))

    double sumCluster =  sumFrictionCoefficient(cluster);
    double sumParticle =  sumFrictionCoefficient(n);

    found_n = false; //TEST
    for(size_t i = 0; i < cluster.size(); i++){
      if(cluster[i] == n)
        found_n = true;
    }
    assert(found_n);

    double P = sumParticle / ( sumCluster );
    //Note, sumParticle is part of sumCluster, since n is part
    //of cluster-vector. This is correct, since then P < 1.
    //or P = 1 if n is not part of a cluster.

    if(1 < P){ //TEST
      cout << "sumParticle" << sumParticle << endl
           << "sumCluster" << sumCluster << endl
           << "Error: P = " << P << endl;
      abort();
    }

    double randomNumb = randomNumber.doub();
    if(randomNumb > P){
      //Don't move the cluster if(P<r):
      //Nothing needs to be done, vacancy-matrix unchanged. If
      //cluster.size()=1, then P=1, and we will always move the particle.
      if(verbose)
        cout << "Cluster to heavy, particle " << n <<" stays \n" << endl;
    }
    else{
      //move whole Cluster.

      if(verbose) cout << "Move the whole cluster, does it collide? ";

      assert(board_.checkOK());

      //We do a trial-and-error move with the whole cluster, since
      //it's the easiest way to check for obstructions/crowders.
      //store the old positions, temporarily:
      vector<int> XoldCluster, YoldCluster, ZoldCluster;
      //TODO: make this use the Particle class beauty...
      for(size_t i = 0; i < cluster.size(); i++){
        int c = cluster[i];
        XoldCluster.push_back(pos_[c].x);
        YoldCluster.push_back(pos_[c].y);
        ZoldCluster.push_back(pos_[c].z);
        //This loop includes particle n.

        //set their old positions as Vacant, perhaps just temporarily:
        board_.setAsVacant(pos_[c]);

        //Move all particles in the cluster:
        moveAndBoundaryCheck(c,direction);
      }

      //Check their new positions against double occupancy.
      //Exclude their old coordinates from the vacancy-matrix.
      //(this is already done by setting them as vacant)

      vector<int> blocking;
      bool occupied = false;
      for(size_t i = 0; i < cluster.size(); i++){
        int c = cluster[i];

        //Now all particles in cluster[] have new positions
        //but their new sites are not yet marked as occupied
        if(board_.isOccupied(pos_[c])){
          int particleLabelOnNewSite = board_.get_label(pos_[c]);
          //store which particles are blocking the cluster;
          blocking.push_back(particleLabelOnNewSite);
          occupied = true;
        }
      }

      //TEST
      for(size_t i = 0; i < blocking.size(); i++ ){
        for(size_t j = 0; j < cluster.size(); j++ ){
          if(cluster[j] == blocking[i])
            throw std::string("Same particle in both vectors!");
        }
      }

      if(verbose){ //is it blocked:
        if(occupied) cout << "YES" << endl << endl;
        else cout << "NO" << endl << endl;
      }

      if(!occupied){
        //no double occupancies, accept move, uppdate vacancy-matrix:

        if(blocking.size() != 0)
          throw std::string("Blocking-vector should have zero length.");

        for(size_t i = 0; i < cluster.size(); i++){
          int c = cluster[i];
          //Now all particles in cluster[] have new positions, including n
          board_.setAsOccupied(pos_[c], c);
        }

        if(verbose){
          cout << "No collisions, moved cluster to:" << endl;
          printPositionToScreen(cluster);
          cout << endl;
        }

        assert(board_.checkOK());

      }
      else{ //Ouch, our cluster is colliding with some shit...
        //NOTE/TODO: vacancy is a bit strange now, with XoldCluster == -1 etc.

        //this is an addition made in July:
        //For fucks sake! I need to restore the vacancy-matrix
        //now, since CountNearestneighbors() and BuildCluster...()
        //depends on it being correct!
        //Reset the "trial-and-error" move (for now), ie.
        //move the cluster back! Pronto!
        for(size_t i = 0; i < cluster.size(); i++){
          int c = cluster[i];
          pos_[c].x = XoldCluster[i];
          pos_[c].y = YoldCluster[i];
          pos_[c].z = ZoldCluster[i];

          board_.setAsOccupied(pos_[c], c);
        }
        assert(board_.checkOK());
        //End of July 2010 addition (but must compensate for this further
        //down when we actually move)
        //NOTE: With the trial-and-error, we got what we wanted, namely
        //the Blocking[]-vector

        assert(!blocking.empty());

        //Construct a vector that contains all path blocking particles.
        buildBlockingCluster(blocking,cluster,direction);

        assert(!blocking.empty());

        if(verbose){
          cout << "The blocking cluster is:" << endl;
          printPositionToScreen(blocking);
          cout << endl;
        }

        double sumBlocking = sumFrictionCoefficient(blocking);

        double P = (double) sumCluster / (sumBlocking + sumCluster);

        double slump = randomNumber.doub();
        if(slump <= P ){
          /*Move Cluster and All the Blocking particles. But beware
            when we move the cluster and update the vacancy matrix we
            must make sure it works for a situation where for instance
            everything gets shifted, around the boundary:
            ---------------------
            | 1 | 2 | 3 | 4 | 5 |
            ---------------------
            if particle 4 is the single particle cluster that pushes
            at the blocking cluster [1 2 3 5], then we must get
            (for X+1 move):
            ---------------------
            | 5 | 1 | 2 | 3 | 4 |
            ---------------------
            and have the vacancy-matrix be correct!
          */

          assert(board_.checkOK());
          assert(!blocking.empty());    //should not be empty

          //Set all old positions as vacant:
          for(size_t i = 0; i < blocking.size(); i++){
            int b = blocking[i];
            board_.setAsVacant(pos_[b]);
          }
          for(size_t i = 0; i < cluster.size(); i++){
            int c = cluster[i];
            board_.setAsVacant(pos_[c]);
          }

          //Move all particles, to get new positions in pos_[].x ...
          for(size_t i = 0; i < blocking.size(); i++){
            int b = blocking[i];
            moveAndBoundaryCheck(b,direction);
          }
          for(size_t i = 0; i < cluster.size(); i++){
            int c = cluster[i];
            moveAndBoundaryCheck(c,direction);
          }

          //Update the vacancy-matrix
          for(size_t i = 0; i < blocking.size(); i++){
            int b = blocking[i];
            board_.setAsOccupied(pos_[b], b);
          }
          for(size_t i = 0; i < cluster.size(); i++){
            int c = cluster[i];
            board_.setAsOccupied(pos_[c], c);
          }

          if(verbose){
            cout << "Move the cluster and the blocking particles." << endl
                 << "Blocking Particles:" << endl;
            printPositionToScreen(blocking);
            cout << "Original Cluster:" << endl;
            printPositionToScreen(cluster);
            cout << endl;
          }
          assert(board_.checkOK());
        }
        else{
          //Move everything back, update vacancy Matrix.

          //SINCE THE JULY ADDITION, WE HAVE ALREADY MOVED BACK THE
          //TRIAL-AND-ERROR MOVE, SO WE NEED NOT DO ANYTHING...

          // for(int i=0; i < Cluster.size(); i++){
          //   int c = Cluster[i];

          //   //Keep in mind: label stored in Cluster[i]
          //   //has its old position stored in X/Y/ZoldCluster[i]
          //   pos_[c].x = XoldCluster[i];
          //   pos_[c].y = YoldCluster[i];
          //   pos_[c].z = ZoldCluster[i];

          //   board_.setAsOccupied(pos_[c], c);
          // }

          if(verbose){
            cout << "The blocking cluster is too heavy, no move is made"
                 << endl << "the cluster is unmoved:" << endl;
            printPositionToScreen(cluster);
            cout << endl;
          }
          assert(board_.checkOK());
        }
        assert(board_.checkOK());
      }//END -- colliding Cluster
      assert(board_.checkOK());
    }//if move whole cluster-END.
    assert(board_.checkOK());
  }//if nearestNeighbor != 0 END
  assert(board_.checkOK());
}


void SuperInteraction::calculateExtent(vector<int>& cluster,int direction,
                                       int& Min, int& Max){
  //Function to calculate the maximum and minimum coordinates of the cluster,
  //perpendicular to the jumping-direction. Only Works for 2D so far...
  //(and 1D? hmm, ingen utbredning vinkelrät mot flyttriktning...)

  //Do I need to consider boundary periodic/fix? I don't think so, since
  //other functions... .... hmmmm

  int dir;
  int clusterSize = cluster.size();

  if( direction == 0 || direction == 1) dir = 0; //Move in X
  if( direction == 2 || direction == 3) dir = 1; //Move in Y
  if( direction == 4 || direction == 5) dir = 2; //Move in Z

  switch(dir){

  case 0: //perpendicular to X
    Max = 1;
    Min = latticeY_;
    for(int i=0; i < clusterSize; i++){
      int m = cluster[i];

      //Max2 = 1;
      //Min2 = latticeZ_;

      if(pos_[m].y < Min ) Min = pos_[m].y;
      if(pos_[m].y > Max ) Max = pos_[m].y;

      //if(pos_[m].z < Min2 ) Min2 = pos_[m].z;
      //if(pos_[m].z > Max2 ) Max2 = pos_[m].z;
    }
    break;


  case 1: //perpendicular to Y
    Max = 1;
    Min = latticeX_;
    for(int i=0; i < clusterSize; i++){
      int m = cluster[i];

      //Max2 = 1;
      //Min2 = latticeZ_;

      if(pos_[m].x < Min ) Min = pos_[m].x;
      if(pos_[m].x > Max ) Max = pos_[m].x;

      //if(pos_[m].z < Min2 ) Min2 = pos_[m].z;
      //if(pos_[m].z > Max2 ) Max2 = pos_[m].z;
    }
    break;

  case 2: //perpendicular to Z
    throw std::string("Move perpendicular to Z, error in 2D!");
    //this code is not yet operational in 3D! --> abort!

    Max = 1;
    Min = latticeX_;
    for(int i=0; i < clusterSize; i++){
      int m = cluster[i];

      //Max2 = 1;
      //Min2 = latticeY_;

      if(pos_[m].x < Min ) Min = pos_[m].x;
      if(pos_[m].x > Max ) Max = pos_[m].x;

      //if(pos_[m].y < Min2 ) Min2 = pos_[m].y;
      //if(pos_[m].y > Max2 ) Max2 = pos_[m].y;
    }
    break;
  }//switch end

  //Check if cluster is wrapped around the periodic boundary
  //-------------------------------------------------------
  bool wrapped = false;
  //move perpendicular to X:
  if(Max == latticeY_ && Min == 1 && dir == 0) wrapped = true;
  //move perpendicular to X: <-- TODO!
  if(Max == latticeX_ && Min == 1 && dir == 1) wrapped = true;
  //if(Max == latticeZ_ && Min == 1 && dir == 3) wrapped = true; //not working, but dir!=3 now...

  if(wrapped){
    //We define a convention:
    //If the cluster is wrapped around the periodic boundary:
    //its borders will be 1 < partOne < MaxX && MinX < partTwo < latticeX_
    //which means that MaxX < MinX, which is how we later will
    //know that the cluster is divided in two.

    //temp set to neg. If still neg in the end -> error!
    int MaxX= -1, MaxY = -1, MinX = -1, MinY = -1;
    bool step = true;

    if(dir == 0){ //If move in X-direction
      //Scan cluster from 1 --> latticeY_ and stop
      // for-loop once step = false,
      for(int i = 1; step && i <= latticeY_; i++){

        step = false; //stop for-loop if no cluster-particles at this site (=i)
        for(int j = 0; j < clusterSize && !step; j++){
          //having "&& !step" in the condition above is just to skip rest of
          //this time-consuming loop over the cluster particles. once we know
          //the cluster has at least one particle here, we don't need to run
          //it more

          if(pos_[cluster[j]].y == i){
            //OK, allow another step upwards in Y
            //to check if the cluster extends this far
            step = true;
            MaxY = i;
          }
        }
      }

      step = true;  //kan nog ta bort denna senare...(när allt funkar)
      //Scan cluster from latticeY_ --> 1
      for(int i = latticeY_; step && i>0; i--){
        step = false;
        for(int j=0; j < clusterSize && !step; j++){
          if(pos_[cluster[j]].y == i){
            //OK, allow another step downwards in Y
            //to check if the cluster extends this far
            step = true;
            MinY = i;
          }
        }
      }
      if(MinY == 1 && MaxY == latticeY_)
        //This is quite possible, since it can be "twisted" like a vortex/swirl
        cout <<"Cluster spans entire lattice!"<<endl;


      Max = MaxY;
      Min = MinY;
    }//"Move in X-direction" -END


    if(dir == 1){ //If move in Y-direction

      //Scan cluster from 1 --> latticeX_
      for(int i = 1; step && i <= latticeX_; i++){
        step = false;
        for(int j = 0; j < clusterSize && !step; j++){
          if(pos_[cluster[j]].x == i){
            //OK, allow another step to the right in X
            //to check if the cluster extends this far
            step = true;
            MaxX = i;
          }
        }
      }

      step = true;  //kan nog ta bort denna senare...(när allt funkar)
      //Scan cluster from latticeX_ --> 1
      for(int i = latticeX_; step && i>0; i--){
        step = false;

        for(int j = 0; j < clusterSize && !step; j++){
          if(pos_[cluster[j]].x == i){
            //OK, allow another step to the left in X
            //to check if the cluster extends this far
            step = true;
            MinX = i;
          }
        }
      }
      if(MinX == 1 && MaxX == latticeX_)
        cout << "Cluster spans entire lattice" << endl;
      //Which doesn't mean there is any error, just some information
      //that's nice to have. Will depend on interaction strength.

      Max = MaxX;
      Min = MinX;
    }//"Move in Y-direction" -END


  }//if "maybe wrapped around periodic"-END
}


void SuperInteraction::buildBlockingCluster(vector<int>& BuildBlockingVector,
                                            vector<int> Cluster, int direction){
  //"Cluster" is the moving cluster that's colliding with BuildBlockingVector
  //"BuildBlockingVector" contains at least one of the blocking particles.
  //We add the rest of the blocking particles to this, and store in unique

  vector<int> unique;


  //TEST to see if Cluster and BuildBlockingVector contains the same particles
  //if so --> check that vacancy matrix hasn't be (incorrectly) changed
  for(size_t i=0; i < BuildBlockingVector.size(); i++ ){
    for(size_t j=0; j < Cluster.size(); j++ ){
      if(Cluster[j] == BuildBlockingVector[i])
        throw std::string("Same particle in both vectors!");
    }
  }

  //First check the min,max coordinates for the Cluster
  //perpendicular to its direction of movement:
  int Min, Max;  //these are given as references in the func. below:
  calculateExtent(Cluster, direction, Min, Max);

  assert(board_.checkOK());

  if(dim_ == 3)
    throw std::string("BuildBlockingCluster does not work in 3D (yet)");

  //if we added a particle to the vector, run again to check its neighbors
  bool CheckAgain;

  unsigned int counter = 0; //debugging variable, see how many times the loop runs
  do{
    counter++;
    CheckAgain = false;

    int limit = BuildBlockingVector.size();
    for(int i=0; i< limit; i++){
      int m = BuildBlockingVector[i];
      board_.checkNeighbours(pos_[m], BuildBlockingVector, isBoundaryFix_);
      //BuildBlockingVector will store more particles each
      //iteration, but only loop over the old particles ( i.e. i < limit)
    }

    for(size_t j = 0; j < BuildBlockingVector.size(); j++){
      bool SingleValued = true;
      int b = BuildBlockingVector[j];

      //Check if it's stored in unique[] already:
      for(size_t k = 0; k < unique.size(); k++){
        if(b == unique[k])
          SingleValued = false;
        //NOTE: Never enters this loop on the first run (unique.size=0)
      }
      //Check that particle s is not part of "Cluster"
      if(SingleValued){ //no need to run it if it's already "false"
        for(size_t i = 0; i < Cluster.size(); i++){
          if(b == Cluster[i])
            SingleValued = false;
        }
      }

      if(SingleValued){

        //Check if Min & Max are in fact wrapped around the boundary!
        bool normal = true;

        //due to definition Max/Min in calculateExtent():
        if(Max < Min) normal = false;

        if(normal && (direction == 0 || direction == 1)){ //Move in X
          if(pos_[b].y <= Max && Min <= pos_[b].y){
            unique.push_back(BuildBlockingVector[j]);
            //New particles added to the cluster, run the loop again,
            //to check the neighbors of this new particle
            CheckAgain = true;
          }
        }
        if(!normal && (direction == 0 || direction == 1)){ //Move in X
          if( (1 <= pos_[b].y && pos_[b].y <=Max ) ||
               ( Min <= pos_[b].y && pos_[b].y <= latticeY_)){
            unique.push_back(BuildBlockingVector[j]);
            //New particles added to the cluster, run the loop again,
            //to check the neighbors of this new particle
            CheckAgain = true;
          }
        }
        if(normal && (direction == 2 || direction == 3)){ //Move in Y
          if(pos_[b].x <= Max && Min <= pos_[b].x){
            unique.push_back(BuildBlockingVector[j]);
            //New particles added to the cluster, run the loop again,
            //to check the neighbors of this new particle
            CheckAgain = true;
          }
        }
        if(!normal && (direction == 2 || direction == 3)){ //Move in Y
          if( (1 <= pos_[b].x && pos_[b].x <=Max ) ||
               ( Min <= pos_[b].x && pos_[b].x <= latticeX_)){
            unique.push_back(BuildBlockingVector[j]);
            //New particles added to the cluster, run the loop again,
            //to check the neighbors of this new particle
            CheckAgain = true;
          }
        }
        if( direction == 4 || direction == 5){ //Move in Z
          throw std::string("I'm in 2D!");
          //throw std::string("I'm in 2D!", __FILE__, __LINE__);
          //this code is not fully functional in 3D....
        }

        /*
        //We have now limited the Blocking cluster on the sides,
        //But not backwards! Exclude any particle that's a part of the Cluster-
        //vector, IF a new particle was added
        if(CheckAgain){
        for(int i = 0; i < Cluster.size(); i++){
        //cout <<endl<<"Cluster.size(): "<<Cluster.size()<<endl;
        //cout <<"i="<<i<<endl;
        //cout << "unique.back: "<<unique.back()<<endl;
        //cout << "unique.size(): "<<unique.size()<<endl<<endl;
        if(unique.back()==Cluster[i]){
        unique.pop_back(); //delete the added particle
        //cout << "unique.size() (after pop_back): "<<unique.size()<<endl;

        //Since we removed the particle, we should not check the neighbors
        //of it in a new run!
        CheckAgain = false;
        }
        }

        //TEST
        if(unique.size()==0){
        cout << "Unique.size() err code 1214: "<<unique.size()<<endl;
        abort();
        }
        }*/

        //TEST
        // if(unique.size()==0){
        //   cout << "unique.size err code 1215: "<<unique.size()<<endl
        //        << "CheckAgain-value: "<< CheckAgain << endl;
        //     abort();
        // }
      }

      //TEST
      // if(unique.size()==0){
      //   cout << "Unique.size() err 1216:"<<unique.size()<<endl;
      //   abort();
      // }


    }

    //TEST
    if(unique.size()==0){
      cout << "Unique.size() err 1217:" << unique.size() << endl;
    }
    //This code works even if BuildBlockingVector-vector contains
    //the same particle multiple times, new and/or old. I think...


    //TEST:
    if(unique.size() > BuildBlockingVector.size() )
      throw std::string("Something wrong in unique[].");

    BuildBlockingVector = unique;

  }while(CheckAgain);

  /*
  //TEST:
  cout << "kollar BuildBlockingVector.size():"<<BuildBlockingVector.size()<<endl;
  for(int j=0; j < BuildBlockingVector.size(); j++){
  int m= BuildBlockingVector[j];
  cout <<"n = "<<m<<" ("<<pos_[m].x<<","<<pos_[m].y<<","<<pos_[m].z<<")"<<endl;
  }

  //TEST
  for(int i=0; i < N; i++){
  if(pos_[i].x==0 || pos_[i].y==0 || pos_[i].z==0){
  cout <<"Utanför brädet!"<<endl;
  abort();
  }
  }
  */

  //TEST to see if Cluster and BuildBlockingVector contains the same particles
  //if so --> check that vacancy matrix hasn't be (incorrectly) changed
  for(size_t i = 0; i < BuildBlockingVector.size(); i++ ){
    for(size_t j = 0; j < Cluster.size(); j++ ){
      if(Cluster[j]==BuildBlockingVector[i])
        throw std::string("Same particle in both vectors (err. location: 2)");
    }
  }

}

void SuperInteraction::buildCluster2(int n, vector<int>& BuildClusterVector,
                                     double Exponent){
  //Build a vector with particles that form the Cluster:
  //-------------------------------------------------
  //each neighbor has a probability given by the Boltzmann factor:
  //    exp(-1.0*Exponent)
  //to belong to the cluster.

  //cout<<"Interaction Strength: "<< Exponent << endl;

  vector<int> unique;
  vector<int> NewParticleAdded;

  double random_number;

  //n included with 100 % probability:
  unique.push_back(n);

  bool CheckAgain;   //true if a new particle was added to the cluster

  do{
    CheckAgain = false;

    int limit = BuildClusterVector.size();
    for(int i = 0; i< limit; i++){
      int m = BuildClusterVector[i];
      board_.checkNeighbours(pos_[m], BuildClusterVector, isBoundaryFix_);
      //BuildClusterVector-vector will store more particles each
      //iteration, but only loop over the first initial set (ie. i < limit)
    }

    for(size_t j = 0; j < BuildClusterVector.size(); j++){
      bool SingleValued = true;

      //Check if it's stored in unique[] already:
      for(size_t k = 0; k < unique.size(); k++){
        if(BuildClusterVector[j] == unique[k])
          SingleValued = false;
      }
      if(SingleValued){
        //Add the new particle, with Boltzmans consent...
        random_number = randomNumber.doub();
        if(random_number > exp(-1.0*Exponent)){
          unique.push_back(BuildClusterVector[j]);
          NewParticleAdded.push_back(BuildClusterVector[j]);
          //New particels added to the cluster, run the loop again
          //for particles in NewParticleAdded to check the neighbors
          // of this new particle
          CheckAgain = true;
        }
      }
    }

    BuildClusterVector = NewParticleAdded;
    NewParticleAdded.clear();

  }while(CheckAgain);

  //Now save the final vector (cluster) and return it as a reference:
  BuildClusterVector = unique;


  if(false){ //deactivated for now
    //TEST:
    cout << "kollar Nearest Neighbour.size(): "
         << BuildClusterVector.size() << endl;
    for(size_t j = 0; j < BuildClusterVector.size(); j++){
      int m= BuildClusterVector[j];
      cout <<"m = " << m << " (" << pos_[m] << ")" << endl;
    }
  }

  //TEST
  for(int i = 0; i < noParticles_; i++){
    if(pos_[i].x == 0 || pos_[i].y == 0 || pos_[i].z == 0)
      throw std::string("Outside the lattice!");
  }

  //TEST---------(To check for double counting)
  bool IsUniqueTrue = true;
  for(size_t i = 0; i < unique.size(); i++){
    for(size_t j = 0; j < unique.size(); j++){
      if(unique[i] == unique[j] && i != j)
        IsUniqueTrue = false;
    }
  }
  if(!IsUniqueTrue)
    throw std::string("Unique is not unique!");
  //-------------
}




void SuperInteraction::moveAndBoundaryCheck(int n, int R){
  //Overriding the default version, since I don't want the
  //vacancy check in this one. (cant do trial and error with
  //entire clusters unless I do this overriding.)

  //Move particle "n" in direction "R", according to the "boundary"
  //rules (fix (1) or periodic (0)), IF the new site is vacant.
  bool boundary = isBoundaryFix_;

  if(0 <= n && n < noParticles_ ){
    if(isTestOn_) cout << "Moving particle " << n << endl;
    // n = particle to move, R = direction to move, such that:
    //      R= 0 = X-right,      1 = X-left
    //         2 = Y-right,      3 = Y-left
    //         4 = Z-up,         5 = Z-down
    // boundary: boundary=0 = periodic, boundary=1 = fix boundary

    //Store in case we want to restore the previous configuration
    Particle old = pos_[n];

    if( R >= 0 && R < 6 ){
      switch(R){
      case 0: pos_[n].x=pos_[n].x + 1;
        if(pos_[n].x == latticeX_ +1)
          pos_[n].x = boundary ? old.x : 1;
        break;
      case 1: pos_[n].x = pos_[n].x - 1; //move particle n to the left
        if(pos_[n].x == 0)
          pos_[n].x =  boundary ? old.x : latticeX_;   //fix : priodic
        break;  //(first coordinate is 1, 0=outside)
      case 2: pos_[n].y = pos_[n].y + 1;
        if(pos_[n].y == latticeY_ +1)
          pos_[n].y = boundary ? old.y : 1;
        break;
      case 3: pos_[n].y = pos_[n].y - 1;
        if(pos_[n].y == 0)
          pos_[n].y = boundary ? old.y : latticeY_;
        break;
      case 4: pos_[n].z = pos_[n].z + 1;
        if(pos_[n].z == latticeZ_+1)
          pos_[n].z = boundary ? old.z : 1;
        break;
      case 5: pos_[n].z = pos_[n].z-1;
        if(pos_[n].z == 0)
          pos_[n].z = boundary ? old.z : latticeZ_;
        break;
      }

      //THIS IS THE DIFFERENCE BETWEEN THE SUPER INTERACTION CODE:
      //(I cant have vacancy check here:)
      //vacancyCheck(n,old);
    }
    else{
      cout << "There is a maximum of " << 2 * dim_
           << " positive directions to move in" << endl;
      abort();
    }
  }
  else{
    cout << "Invalid particle label, n=" << n
         << ", must fulfill 0 <= n < N" << endl;
    abort();
  }
}



//Override the lattice-class version:
void SuperInteraction::move(){

  double tau;      //waiting time
  int i = 0;       //index of samplingTime-vector
  timeSum_ = 0;

  while(timeSum_ < samplingTime_[noSamplingTimes_-1]){
    //"-1" since we start on 0
    //The time-sampling was completely rewritten in late August 2010
    //to resolve the spaghetti that was the previous version. This
    //follows the Gillespie_exclusion2.cpp-implementation closely, to
    //make it easier on those already familiar with that source code.

    tau = computeWaitingTime();

    //Save displacement if next time-step is beyond next sampling time
    while(timeSum_ <= samplingTime_[i] &&
          samplingTime_[i] < timeSum_ + tau  && i < noSamplingTimes_){

      //save displacement (from previous step)
      dx_[i] = pos_[0].x - pos_0_.x;
      dy_[i] = pos_[0].y - pos_0_.y;
      dz_[i] = pos_[0].z - pos_0_.z;
      dr_[i] = sqrt( pow(dx_[i],2) + pow(dy_[i],2) + pow(dz_[i],2) );

      i++;
    }

    /*
    //TEST (this is similar to what my spaghetti code used:)
    if(timeSum_ <= samplingTime_[i] &&  samplingTime_[i] < timeSum_ + tau &&
    i < noSamplingTimes_){

    //save displacement (from previous step)
    dx[i] = pos_[0].x - pos_0.x;
    dy[i] = pos_[0].y - pos_0.y;
    dz[i] = pos_[0].z - pos_0.z;
    dr[i] = sqrt( pow(dx[i],2) + pow(dy[i],2) + pow(dz[i],2) );

    i++;
    }
    */


    //TODO remove! -------------------------------------------
    //pick a particle & direction based on jump-rate, and move.
    if(partialSum_.empty())
      throw std::string("Partial-sum vector has not been initiated");

    double r2;                          //store random number here

    int mu_guess;                       // must be integer. (index of vector)
    int mu_left = 0;
    int mu_right = dim_*2*noParticles_;
    double p_left = (double) partialSum_.front();  //value in first element ( =0 )
    double p_right = (double) partialSum_.back();  //value in last element ( =dim_*2*N)

    do{
      r2 = randomNumber.doub();
    }while( r2 == 1 || r2 == 0 );

    double p_rand = (double) r2*partialSum_.back();
    bool loopAgain = true;

    //NOTE: when comparing a long double against a double, C++ will fill
    //in the blanks with zeroes. Every 10^15 run of move() will lead to
    //a situation where we actually have p_rand==p[mu_guess], and
    //changing the interval below in the if statements,
    //from "..<= .. <.." to "..< .. <=.." will result in an infinite
    //loop every 10^15 turns since "double" stores 52 bit, 2^52 = 10^15.
    //This bug has now been fixed!

    //finds the mu
    do{
      mu_guess = (int) ( (p_rand - p_left) * ( mu_right - mu_left )
                         /(p_right - p_left) + mu_left );

      if((double) partialSum_[mu_guess] <= p_rand && p_rand <
          (double) partialSum_[mu_guess+1])
        loopAgain = false;

      else{
        if(p_rand < (double) partialSum_[mu_guess]){
          mu_right = mu_guess;
          p_right = partialSum_[mu_guess];
        }
        else{//if(r2*partialSum_.back() > partialSum_[mu_guess+1]){
          mu_left = mu_guess + 1;
          p_left = partialSum_[mu_guess+1];
        }
      }
    }while(loopAgain);

    int mu = mu_guess;
    int n,r;

    //transforms the index mu to which particle to move
    //(index n, by reference) and in which direction,by returning
    //an integer 0<= r <= 5 (in 3D)
    convertMuToParticle(mu,n,r);

    //use this one instead of Lattice::moveAndBoundaryCheck():
    superInteractionCode(n,r);

    timeSum_ = timeSum_ + tau;

  }
}
