// Guy Shechter
// June 2004

// Variables used for time debugging.
TV tv1,tv2; 
TZ tz;

//
// Call this function to build the k-d Tree.
//
Tree *build_kdtree(double *reference, int N, int D, int *index, 
		   int L,int offset) {

  Tree *treeptr;
  if ( (treeptr = (Tree *) malloc(sizeof(Tree))) == NULL )
    mexPrintf("Error allocating memory for kD-Tree\n");
  
  treeptr->rootptr = build_kdtree_core(reference, N, D, index, L, offset);
  treeptr->dims = D;

  return treeptr;
}

//
// This is the core function for building the k-D tree for raw input data.
//
Node *build_kdtree_core(double *reference, int N, int D, int *index, 
			int L,int offset) {
  int i,median;
  Node *nodeptr;
  
#ifdef DEBUG_BUILD_TREE
  mexPrintf("Entering build_kdtree \n");
#endif
  
  if ( (nodeptr = (Node *) malloc(sizeof(Node))) == NULL )
    mexPrintf("Error allocating memory for kD-Tree\n");

  if ( (nodeptr->pt = (double*) malloc(sizeof(double)*D)) == NULL )
    mexPrintf("Error allocating memory for kD-Tree\n");
  
  if (L==1) {
    for (i=0; i < D; i++)
      nodeptr->pt[i] = EVAL_INDEX(index[0],i,N);
    nodeptr->index = (unsigned int) index[0];
    
#ifdef DEBUG_BUILD_TREE
    mexPrintf("Created leaf     node, (dim=%d) (value: ", offset);
    for (i=0; i < D; i++)
      mexPrintf(" %f",nodeptr->pt[i]);
    mexPrintf(" )\n");
#endif
    
    nodeptr->orientation = offset;
    nodeptr->left=NULL;
    nodeptr->right=NULL;
    return nodeptr;
  }
  
  quicksort(index,0,L-1,reference,offset,N);
  median = (int)((L)/2.0);
  
  for (i=0; i < D; i++)
    nodeptr->pt[i] = EVAL_INDEX(index[median],i,N);
  nodeptr->index = (unsigned int) index[median];
  nodeptr->orientation = offset;
  nodeptr->left=NULL;
  nodeptr->right=NULL;
  
#ifdef DEBUG_BUILD_TREE
  mexPrintf("Created internal node, (dim=%d) (value: ", offset);
  for (i=0; i < D; i++)
    mexPrintf(" %f",nodeptr->pt[i]);
  mexPrintf(" )\n");
#endif
  
  nodeptr->left = build_kdtree_core(reference,N,D,index,median,(offset+1)%D);
  
  if (L-median-1)
    nodeptr->right = build_kdtree_core(reference,N,D,index+median+1,
				       L-median-1,(offset+1)%D);
  
  return nodeptr;
}


// 
//   The following two functions are for the quicksort algorithm 
//
int partition(int *a, int p, int r, double *reference, int offset, int D){
  int i,j;
  double x,tmp;
  
  x= EVAL_INDEX(a[p],offset,D); i=p-1; j=r+1;
  while (1) {
    for (j--; EVAL_INDEX(a[j],offset,D) > x; j--);
    for (i++; EVAL_INDEX(a[i],offset,D) <x; i++);
    if (i<j){
      INT_SWAP(a[j],a[i]);
    }
    else return j;
  }
}

void quicksort(int *ra, int p, int r, double *reference, int offset,int D){
  int q;
  
  if (p < r) {
    q = partition(ra,p,r,reference,offset,D);
    quicksort(ra,p,q,reference,offset,D);
    quicksort(ra,q+1,r,reference,offset,D);
  }
}



//
// Deallocate memory used by the k-D tree
//
void free_tree(Node *pVertex) {
  
  if (pVertex == NULL)
    return;
  if (pVertex->left)  free_tree(pVertex->left);
  if (pVertex->right) free_tree(pVertex->right);
  
  free(pVertex->pt);
  free(pVertex);
  return;
}


//
// Compute the norm distance between two points of dimension Dim
//
double calcdistance(double *pt1, double *pt2, int Dim){
  int j;
  double d=0;
  for (j=0; j < Dim; j++) 
    d+= (pt1[j]-pt2[j])*(pt1[j]-pt2[j]);
  return (sqrt(d));
}


Node *pointLocation(Node *v, double *pt,int D){
  
  //if (!v) v=root;

  if ( ( v->left == NULL) && (v->right == NULL) ) 
    return v;

  if (    pt[v->orientation] < v->pt[ v->orientation] )
    return ( (v->left) ? pointLocation(v->left,pt,D) : v);
  else if (    pt[v->orientation] > v->pt[v->orientation]  ) 
    return ( (v->right) ? pointLocation(v->right,pt,D): v);
  else {
    //What we have is pt[v->orientation] == v->pt[v->orientation] 
    Node *vl = ( (v->left)? pointLocation(v->left,pt,D) : v);
    Node *vr = ( (v->right)? pointLocation(v->right,pt,D) : v);
    if ( calcdistance(vl->pt,pt,D) < calcdistance(vr->pt,pt,D) )
      return vl;
    else
      return vr;
  }
}


Node * rangeQuery(Node *v, double distance, double *pt, int D){

  Node *current=NULL, *tmp;
  int j;

  
  if ( ( (v->pt[v->orientation] - distance) <=  pt[v->orientation]) &&
       ( (v->pt[v->orientation] + distance) >=  pt[v->orientation]) ){
    if (calcdistance(pt,v->pt,D) <= distance){
      current=(Node*) malloc(sizeof(Node));
      current->pt=(double*) malloc(sizeof(double)*D);
      current->right=NULL;
      for (j=0; j < D; j++) current->pt[j]=v->pt[j];
      current->index=v->index;
    }
  }

  
  if (  (!(v->left)) && (!(v->right)) ) 
    return current;
  
  tmp=current;
  while (tmp && tmp->right) 
    tmp=tmp->right;
  
  if (v->left){
    if ( v->pt[v->orientation] >= (pt[v->orientation] - distance ) ){
      if (tmp)
	tmp->right=rangeQuery(v->left, distance, pt,D);
      else
	current=rangeQuery(v->left, distance, pt, D);
    }
  }
  
  tmp=current;
  while (tmp &&tmp->right) 
    tmp=tmp->right;
  
  if (v->right){
    if (   v->pt[v->orientation]  <= (pt[v->orientation] + distance ) ){
      if (tmp)
	tmp->right=rangeQuery(v->right, distance, pt, D);
      else
	current=   rangeQuery(v->right, distance, pt, D);
    }
  }    

  return current;
}

//
// The top-level function for running a query on the k-D tree.
//
void run_queries( Node *pVertex, double *model, int M, int D, 
		  double *closest_pt, double *distance, short ReturnType) {
  
  int i,j;
  double min_distance, *pt;
  Node *LL, *cur, *leaf, *tmp;
  
  pt= (double*)malloc(sizeof(double)*D);
  
  for (i=0; i < M; i++) {
    
#ifdef DEBUG_RUN_QUERIES
    mexPrintf("Running Query (%d/%d) (value: ",
	      i+1, M);
    for (j=0; j < D ; j++)
      mexPrintf(" %f", model[ M*j+i]);
    mexPrintf(" )\n");
#endif
    
    for (j=0; j < D; j++) 
      pt[j]=model[M*j+i];
    
    leaf=pointLocation(pVertex,pt,D);
    min_distance=calcdistance(leaf->pt, pt,D )+0.001;

    LL=rangeQuery(pVertex,min_distance,pt,D);

    if (!LL) {
      if (ReturnType == RETURN_INDEX) 
	closest_pt[i] = -1;
      else{
	for (j=0; j< D; j++)
	  closest_pt[j*M+i]=-1;
      }
      mexPrintf("Null LL\n");
    } 
    else {
      distance[i]=calcdistance(LL->pt, pt,D);
      if (ReturnType == RETURN_INDEX) 
	closest_pt[i] = LL->index;
      else {
	for (j=0; j < D; j++) 
	  closest_pt[j*M+i] = LL->pt[j];
      }
      cur=LL;
      while (cur){
	if ( calcdistance(cur->pt, pt,D) <= distance[i] ) {
	  if (ReturnType == RETURN_INDEX) 
	    closest_pt[i] = cur->index;
	  else {
	    for (j=0; j < D; j++) 
	      closest_pt[j*M+i] = cur->pt[j];
	  }
	  distance[i]=calcdistance(cur->pt, pt,D);
	}
	tmp=cur;
	cur=cur->right;
	free(tmp->pt);
	free(tmp);
      }
    }

#ifdef DEBUG_RUN_QUERIES
    mexPrintf("Distance to closest point is %f\n",distance[i]);
#endif

  }
  free(pt);
}



// 
// Outputs the k-D tree while performing a depth first traversal.
//
void display_tree(Node *nodeptr,int D) {
  int i;
  
  mexPrintf("Node: (dim = %d) (idx= %u) (value = ", nodeptr->orientation,
	    nodeptr->index);
  for (i=0; i <D; i++)
    mexPrintf("%f\t",nodeptr->pt[i]);
  mexPrintf(")\n");

  if (nodeptr->left == NULL) 
    mexPrintf("Left is (null)\n");
  else{
    mexPrintf("Going left\n");
    display_tree(nodeptr->left,D);
  }
  
  if (nodeptr->right == NULL) 
    mexPrintf("Right is (null)\n");
  else{
    mexPrintf("Going right\n");
    display_tree(nodeptr->right,D);
  }
}









//
// The top-level function for running a range search on the k-D tree.
//
void run_range_search( Node *pVertex, double *model, int M, int D, 
		       double distlim, double **pts_in_range, 
		       unsigned int *L, unsigned int **indices){
  
  int i,j, count;
  double min_distance;
  double *pt;
  Node *LL, *cur, *leaf, *tmp;

  pt= (double*)malloc(sizeof(double)*D);

  for (i=0; i < M; i++) {

#ifdef DEBUG_RUN_RANGE_SEARCH
    mexPrintf("Running Search (%d/%d) (value: ",
	      i+1, M);
    for (j=0; j < D ; j++)
      mexPrintf(" %f", model[ M*j+i]);
    mexPrintf(" )\n");
#endif

    for (j=0; j < D ; j++)
      pt[j]=model[M*j+i];
    LL=rangeQuery(pVertex,distlim,pt,D);

    if (!LL) {
      *L=0;
      (*pts_in_range)=NULL;
#ifdef DEBUG_RUN_RANGE_SEARCH
      mexPrintf("Null LL\n");
#endif
    } 
    else {
      cur=LL;
      count=0;
      while (cur){
	cur=cur->right; count++;
      }
#ifdef DEBUG_RUN_RANGE_SEARCH
      mexPrintf("Found %d points\n",count);
#endif
      *L=count;
      (*indices) = (unsigned int *) malloc (sizeof(unsigned int)*count);
      (*pts_in_range) = (double *) malloc (sizeof(double) *count*D);


      cur=LL;
      count=0;
      while(cur){
	(*indices)[count] = cur->index;
	for (j=0; j < D; j++) {
	  (*pts_in_range)[j*(*L)+count] = cur->pt[j];
	}
	tmp=cur;
	cur=cur->right;
	free(tmp->pt);
	free(tmp);
	count++;
      }
    }
  }
  free(pt);
}

