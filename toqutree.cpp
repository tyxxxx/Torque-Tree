
/**
 *
 * toqutree (pa3)
 * significant modification of a quadtree .
 * toqutree.cpp
 * This file will be used for grading.
 *
 */

#include "toqutree.h"

toqutree::Node::Node(pair<int,int> ctr, int dim, HSLAPixel a)
		:center(ctr),dimension(dim),avg(a),NW(NULL),NE(NULL),SE(NULL),SW(NULL)
{}

toqutree::~toqutree(){
	clear(root);
}

toqutree::toqutree(const toqutree & other) {
	root = copy(other.root);
}


toqutree & toqutree::operator=(const toqutree & rhs){
	if (this != &rhs) {
		clear(root);
		root = copy(rhs.root);
	}
	return *this;
}

toqutree::toqutree(PNG & imIn, int k){
	int ul_x = (int) (imIn.width() - pow(2,k))/2;
	int ul_y = (int) (imIn.height() - pow(2,k))/2;
	auto length = (unsigned int) pow(2,k);
	PNG im = PNG(length,length);
	for (unsigned int x = 0; x < length; x++) {
		for (unsigned int y = 0; y < length; y++) {
			*(im.getPixel(x,y)) = *imIn.getPixel(ul_x+x, ul_y+y);
		}
	}
	root = buildTree(&im,k);

/* This constructor grabs the 2^k x 2^k sub-image centered */
/* in imIn and uses it to build a quadtree. It may assume  */
/* that imIn is large enough to contain an image of that size. */

/* your code here */
}

int toqutree::size() {
	return Num(root);

}

int toqutree::Num(Node*curr){
	if(curr){
		return 1+Num(curr->NW)+Num(curr->NE)+Num(curr->SE)+Num(curr->SW);
	} else {
		return 0;
	}

}


toqutree::Node * toqutree::buildTree(PNG * im, int k) {

	if (k >= 0){
		int ul_x = (int) (im->width() - pow(2,k))/2;
		int ul_y = (int) (im->height() - pow(2,k))/2;
		int lr_x = ul_x + (int) pow(2,k)-1;
		int lr_y = ul_y + (int) pow(2,k)-1;
		pair<int,int> ul = make_pair(ul_x,ul_y);
		pair<int,int> lr = make_pair (lr_x,lr_y);
		stats* s =  new stats(*im);
		HSLAPixel avg = s->getAvg(ul,lr);
		pair<int ,int > center = findSplitPoint(im, k);
		Node* root = new Node(center,k,avg);
    int length = (int) pow(2,k);
		int quadrant_length = (int) pow(2,k)/2;
		PNG* NW_im = TreeHelper(im, (center.first+quadrant_length)%length, (center.second+quadrant_length)%length);
		PNG* NE_im = TreeHelper(im, center.first, (center.second+quadrant_length)%length);
		PNG* SE_im = TreeHelper(im, center.first, center.second);
		PNG* SW_im = TreeHelper(im, (center.first+quadrant_length)%length, center.second);
		root->NW = buildTree(NW_im, k-1);
		root->NE = buildTree(NE_im, k-1);
		root->SE = buildTree(SE_im, k-1);
		root->SW = buildTree(SW_im, k-1);
		delete s;
		delete NW_im;
		delete NE_im;
		delete SE_im;
		delete SW_im;
		return root;

	} else{
		return nullptr;
	}

/* your code here */


}


PNG* toqutree::TreeHelper(PNG *im, int x, int y){
	unsigned int Flength = im->width();
	unsigned int length = (im->width())/2;
	PNG* AnotherIm = new PNG(length,length);
	for (unsigned int i = 0; i < length; i++) {
		for (unsigned int j = 0; j < length; j++) {
			*(AnotherIm->getPixel(i,j)) = *(im->getPixel((x+i)%Flength, (y+j)%Flength));
		}
	}
	return AnotherIm;
}

pair<int,int> toqutree::findSplitPoint(PNG* im,int k) {

	if (k == 0) {
		return make_pair(0, 0);
	}
	if (k == 1) {
		return make_pair(1, 1);
	}
	
	int width = im->width();
	int height = im->height();
	int length = pow(2, k - 1);
	double minEntropy = INT32_MAX;
	double entropySum = 0;
	double entropySE = 0;
	double entropyNE = 0;
	double entropySW = 0;
	double entropyNW = 0;
	stats* s = new stats(*im);
	pair<int, int> center;
	// = make_pair(0, 0);
	//3 * pow(2, k) / 4 - 1 
	//pow(2, k) + length - 1
	for (int i = pow(2, k) / 4; i <= pow(2, k)/4 + length - 1; i++) {
		for (int j = pow(2, k) / 4; j <= pow(2, k)/4 + length - 1; j++) {
			int ul_x = i;
			int ul_y = j;
			int lr_x = ul_x + length - 1;
			int lr_y = ul_y + length - 1;
			//when lr_y + 1 < height
			if (lr_y + 1 < height) {
				if (lr_x + 1 < width) {
					entropySE = s->entropy(make_pair(ul_x, ul_y), make_pair(lr_x, lr_y));
					entropyNE = findEntropy2(make_pair(ul_x, lr_y + 1), make_pair(lr_x, height - 1), make_pair(ul_x, 0), make_pair(lr_x, ul_y - 1), k - 1, s);
					entropySW = findEntropy2(make_pair(lr_x + 1, ul_y), make_pair(width - 1, lr_y), make_pair(0, ul_y), make_pair(ul_x - 1, lr_y), k - 1, s);
					entropyNW = findEntropy4(make_pair(lr_x + 1, lr_y + 1), make_pair(width - 1, height - 1), make_pair(0, 0), make_pair(ul_x - 1, ul_y - 1),
							make_pair(0, lr_y + 1), make_pair(ul_x - 1, height - 1), make_pair(lr_x + 1, 0), make_pair(width - 1, ul_y - 1), k - 1, s);
				}
				else if (lr_x + 1 == width) {
					entropySE = s->entropy(make_pair(ul_x, ul_y), make_pair(lr_x, lr_y));
					entropyNE = findEntropy2(make_pair(ul_x, lr_y + 1), make_pair(lr_x, height - 1), make_pair(ul_x, 0), make_pair(lr_x, ul_y - 1), k - 1, s);
					entropySW = s->entropy(make_pair(0, ul_y), make_pair(ul_x - 1, lr_y));
					entropyNW = findEntropy2(make_pair(0, lr_y + 1), make_pair(ul_x - 1, height - 1), make_pair(0, 0), make_pair(ul_x - 1, ul_y - 1), k - 1, s);
				}
				else if (lr_x + 1 > width) {
					lr_x = lr_x % width;
					entropySE = findEntropy2(make_pair(ul_x, ul_y), make_pair(width - 1, lr_y), make_pair(0, ul_y), make_pair(lr_x, lr_y), k - 1, s);
					entropyNE = findEntropy4(make_pair(ul_x, lr_y + 1), make_pair(width - 1, height - 1), make_pair(0, 0), make_pair(lr_x, ul_y - 1),
							make_pair(0, lr_y + 1), make_pair(lr_x, height - 1), make_pair(ul_x, 0), make_pair(width - 1, ul_y - 1), k - 1, s);
					entropySW = s->entropy(make_pair(lr_x + 1, ul_y), make_pair(ul_x - 1, lr_y));
					entropyNW = findEntropy2(make_pair(lr_x + 1, lr_y + 1), make_pair(ul_x - 1, height - 1), make_pair(lr_x + 1, 0), make_pair(ul_x - 1, ul_y - 1), k - 1, s);
				}
			}
			// when lr_y + 1 == height
			else if (lr_y + 1 == height) {
				if (lr_x + 1 < width) {
					entropySE = s->entropy(make_pair(ul_x, ul_y), make_pair(lr_x, lr_y));
					entropyNE = s->entropy(make_pair(ul_x, 0), make_pair(lr_x, ul_y - 1));
					entropySW = findEntropy2(make_pair(lr_x + 1, ul_y), make_pair(width - 1, lr_y), make_pair(0, ul_y), make_pair(ul_x - 1, lr_y), k - 1, s);
					entropyNW = findEntropy2(make_pair(lr_x + 1, 0), make_pair(width - 1, ul_y - 1), make_pair(0, 0), make_pair(ul_x - 1, ul_y - 1), k - 1, s);
				}
				else if (lr_x + 1 == width) {
					entropySE = s->entropy(make_pair(ul_x, ul_y), make_pair(lr_x, lr_y));
					entropyNE = s->entropy(make_pair(ul_x, 0), make_pair(lr_x, ul_y - 1));
					entropySW = s->entropy(make_pair(0, ul_y), make_pair(ul_x - 1, lr_y));
					entropyNW = s->entropy(make_pair(0, 0), make_pair(ul_x - 1, ul_y - 1));
				}
				else if (lr_x + 1 > width) {
					lr_x = lr_x % width;
					entropySE = findEntropy2(make_pair(ul_x, ul_y), make_pair(width - 1, height - 1), make_pair(0, ul_y), make_pair(lr_x, height - 1),k-1,s);
					entropyNE = findEntropy2(make_pair(ul_x, 0), make_pair(width - 1, ul_y - 1), make_pair(0, 0), make_pair(lr_x, ul_y - 1),k-1,s);
					entropySW = s->entropy(make_pair(lr_x + 1, ul_y), make_pair(ul_x - 1, height - 1));
					entropyNW = s->entropy(make_pair(lr_x + 1, 0), make_pair(ul_x - 1, ul_y - 1));
				}
			}
	        // when lr_y + 1 > height
			else if (lr_y + 1 > height) {
				lr_y = lr_y % height;
				if (lr_x + 1 < width) {
					entropySE = findEntropy2(make_pair(ul_x, ul_y), make_pair(lr_x, height - 1), make_pair(ul_x, 0), make_pair(lr_x, lr_y),k-1,s);
					entropyNE = s->entropy(make_pair(ul_x, lr_y + 1), make_pair(lr_x, ul_y - 1));
					entropySW = findEntropy4(make_pair(lr_x + 1, ul_y), make_pair(width - 1, height - 1), make_pair(0, 0), make_pair(ul_x - 1, lr_y),
							make_pair(0, ul_y), make_pair(ul_x - 1, height - 1), make_pair(lr_x + 1, 0), make_pair(width - 1, lr_y),k-1,s);
					entropyNW = findEntropy2(make_pair(lr_x + 1, lr_y + 1), make_pair(width - 1, ul_y - 1), make_pair(0, lr_y + 1), make_pair(ul_x - 1, ul_y - 1),k-1,s);
				}
				else if (lr_x + 1 == width) {
					entropySE = findEntropy2(make_pair(ul_x, ul_y), make_pair(width - 1, height - 1), make_pair(ul_x, 0), make_pair(lr_x, lr_y), k - 1, s);
					entropyNE = s->entropy(make_pair(ul_x, lr_y + 1), make_pair(lr_x, ul_y - 1));
					entropySW = findEntropy2(make_pair(0, ul_y), make_pair(ul_x - 1, height - 1), make_pair(0, 0), make_pair(ul_x - 1, lr_y),k-1,s);
					entropyNW = s->entropy(make_pair(0, lr_y + 1), make_pair(ul_x - 1, ul_y - 1));
				}
				else if (lr_x + 1 > width) {
					lr_x = lr_x % width;
					entropySE = findEntropy4(make_pair(ul_x, ul_y), make_pair(width - 1, height - 1), make_pair(0, 0), make_pair(lr_x, lr_y),
							make_pair(0, ul_y), make_pair(lr_x, height - 1), make_pair(ul_x, 0), make_pair(width - 1, lr_y),k-1,s);
					entropyNE = findEntropy2(make_pair(ul_x, lr_y + 1), make_pair(width - 1, ul_y - 1), make_pair(0, lr_y + 1), make_pair(lr_x, ul_y - 1),k-1,s);
					entropySW = findEntropy2(make_pair(lr_x + 1, ul_y), make_pair(ul_x - 1, height - 1), make_pair(lr_x + 1, 0), make_pair(ul_x - 1, lr_y),k-1,s);
					entropyNW = s->entropy(make_pair(lr_x + 1, lr_y + 1), make_pair(ul_x - 1, ul_y - 1));
				}
			}
			

	//		entropySE = findEntropy(im, make_pair(i, j), k - 1);
	//		entropyNE = findEntropy(im, make_pair(i, (j + length) % size), k - 1);
	//		entropySW = findEntropy(im, make_pair((i + length) % size, j), k - 1);
	//		entropyNW = findEntropy(im, make_pair((i + length) % size, (j + length) % size), k - 1);


			entropySum = entropySE + entropyNE + entropySW + entropyNW;
			if (entropySum < minEntropy) {
				minEntropy = entropySum;
				center = make_pair(i, j);
			}
		}
	}
	delete s;
	return center;

}

//double toqutree::findEntropy(PNG* im, pair<int,int> ul, int k) {
//	int parentLen = pow(2, k + 1);
//	int length = pow(2, k);
//	PNG* subImage = new PNG(length, length); 
//	for (int i = 0; i < length; i++) {
//		for (int j = 0; j < length; j++) {
//			*(subImage->getPixel(i, j)) = *(im->getPixel((ul.first + i) % parentLen, (ul.second + j) % parentLen));
//		}
//	}
//	stats* s = new stats(*subImage);
//	double ret = s->entropy(make_pair(0, 0), make_pair(length - 1, length - 1));
//	delete subImage;
//	delete s;
//	return ret;

//	vector<int> hist;

//}

double toqutree::findEntropy2(pair<int, int> ul1, pair<int, int>lr1, pair<int, int> ul2, pair<int, int>lr2, int k, stats* s) {
	int length = pow(2, k);
	vector<int> hist;
	vector<int> b1Hist;
	vector<int> b2Hist;
	b1Hist = s->buildHist(ul1, lr1);
	b2Hist = s->buildHist(ul2, lr2);
	for (int i = 0; i < 36; i++)
		hist.push_back(b1Hist[i] + b2Hist[i]);
//		b1Hist[i] += b2Hist[i];
	double ret = s->entropy(hist, length*length);
	return ret;
//	return s->entropy(b1Hist, length*length);
}

double toqutree::findEntropy4(pair<int, int> ul1, pair<int, int>lr1, pair<int, int> ul2, pair<int, int>lr2, 
	pair<int, int> ul3, pair<int, int>lr3, pair<int, int> ul4, pair<int, int>lr4, int k, stats* s) {
	int length = pow(2, k);
	vector<int> hist;
	vector<int> b1Hist;
	vector<int> b2Hist;
	vector<int> b3Hist;
	vector<int> b4Hist;
	b1Hist = s->buildHist(ul1, lr1);
	b2Hist = s->buildHist(ul2, lr2);
	b3Hist = s->buildHist(ul3, lr3);
	b4Hist = s->buildHist(ul4, lr4);
	for (int i = 0; i < 36; i++)
		hist.push_back(b1Hist[i] + b2Hist[i] + b3Hist[i] + b4Hist[i]);
//		b1Hist[i] += b2Hist[i] + b3Hist[i] + b4Hist[i];
	double ret = s->entropy(hist, length*length);
	return ret;
	
	
//	return s->entropy(b1Hist, length*length);
}




PNG toqutree::render(){
	int dim = root->dimension;
	auto length = pow(2,dim);
	PNG png = PNG(length,length);
	render(root, &png);
	return png;


/* your code here */

}

void toqutree::render(Node *curr, PNG *png){
	int width =(int) pow(2,curr->dimension);

	if(curr->SE==NULL){
		for(int i=0;i<width;i++){
			for(int j=0;j<width;j++){
				*(png->getPixel(i,j))=curr->avg;
			}
		}
		return;

	}else{
		PNG *newIm= new PNG(width/2,width/2);
		int splitX=curr->center.first;
		int splitY=curr->center.second;
		render(curr->SE,newIm);
		for(int i=0;i<width/2;i++){
			for(int j=0;j<width/2;j++){
				*(png->getPixel((splitX+i)%width,(splitY+j)%width))=*(newIm->getPixel(i,j));
			}
		}
		render(curr->SW,newIm);
		for(int i=0;i<width/2;i++){
			for(int j=0;j<width/2;j++){
				*(png->getPixel((splitX+width/2+i)%width,(splitY+j)%width))=*(newIm->getPixel(i,j));
			}
		}
		render(curr->NE,newIm);
		for(int i=0;i<width/2;i++){
			for(int j=0;j<width/2;j++){
				*(png->getPixel((splitX+i)%width,(splitY+width/2+j)%width))=*(newIm->getPixel(i,j));
			}
		}
		render(curr->NW,newIm);
		for(int i=0;i<width/2;i++){
			for(int j=0;j<width/2;j++){
				*(png->getPixel((splitX+width/2+i)%width,(splitY+width/2+j)%width))=*(newIm->getPixel(i,j));
			}
		}
		delete newIm;
	}
}


void toqutree::prune(double tol){
	pruneHelper1(root,tol);
}

void toqutree::pruneHelper1(toqutree::Node *node, double tol) {
	if (node -> SW == NULL) {
		return;
	}
	if (pruneHelper2(node, node, tol)){
		clear(node->SW);
		clear(node->SE);
		clear(node->NE);
		clear(node->NW);
	} else {
		pruneHelper1(node->SW, tol);
		pruneHelper1(node->SE, tol);
		pruneHelper1(node->NE, tol);
		pruneHelper1(node->NW, tol);
	}
}

bool toqutree::pruneHelper2(Node* curr, Node* curr1, double tol){
	if (curr->SW == NULL){
		return curr->avg.dist(curr1->avg) <= tol;
	} else {
		return pruneHelper2(curr->SW, curr1, tol) && pruneHelper2(curr->SE, curr1, tol) && pruneHelper2(curr->NE, curr1, tol) && pruneHelper2(curr->NW, curr1, tol);
	}
}


void toqutree::clear(Node * & curr){
	if (curr==NULL) {
		return;
	}
	clear(curr->SE);
	clear(curr->NE);
	clear(curr->NW);
	clear(curr->SW);
	delete curr;
	curr = NULL;
/* your code here */
}

/* done */

toqutree::Node * toqutree::copy(const Node * other) {
	if (other == NULL) return NULL;
	Node* temp = new Node(other->center,other->dimension,other->avg);
	temp->SW = copy(other->SW);
	temp->NW = copy(other->NW);
	temp->NE = copy(other->NE);
	temp->SE = copy(other->SE);
	return temp;
/* your code here */
}
