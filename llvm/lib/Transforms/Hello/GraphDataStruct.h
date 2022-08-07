#include <iostream>
using std::cout;
using std::ostream;

#include <string>
using std::string;

#include <vector>
using std::vector;
#include <algorithm>

#include <memory>
using std::unique_ptr;

#include <stack>
using std::stack;

class Vertex {

	public:
		string label;
		vector<Vertex*> outgoingEdges;

	public:

		Vertex(const string & labelIn)
			: label(labelIn) {
			}

		void addEdgeTo(Vertex * to) {
			outgoingEdges.push_back(to);
		}

		void write(llvm::raw_fd_ostream & o) const {

			o << label ;//<< " -> [";
			//for (auto v : outgoingEdges) {
			//	o << " " << v->label;
			//}
			//o << " ]";

		}
};

llvm::raw_fd_ostream & operator<<(llvm::raw_fd_ostream & o, const Vertex & toPrint) {
	toPrint.write(o);
	return o;
}


class DiGraph {

	public:
		vector<unique_ptr<Vertex> > vertices;

	public:

		void addVertex(const string & labelIn) {
			//addVertex.push_back(Vertex(labelIn));
      for(auto const & v : vertices){
        if(v->label == labelIn)
          return;
      }

			vertices.emplace_back(new Vertex(labelIn));
		}

		Vertex* operator[](const int i) {
			return vertices[i].get();
		}

		Vertex * findVertexByName(string name){
			for(int i = 0; i < vertices.size(); i++){
				if(vertices[i] -> label == name)
					return &(*vertices[i]);
			}
		}
		int findVertexIndex(string name){
			for(int i = 0; i < vertices.size(); i++){
				if(vertices[i] -> label == name)
					return i;
			}
		}
		void addEdgeTo(Vertex * from, Vertex * to){
			from->addEdgeTo(to);
		}
		void addEdgeTo(string from, string to){
			findVertexByName(from) -> addEdgeTo(findVertexByName(to));
		}


//function definition
bool isPathDFS(std::unique_ptr<Vertex>  v, std::unique_ptr<Vertex> w){
  return isPathDFS(v->label, w->label);
}
bool isPathDFS(string source, string end){
    //DFS implementation 
    bool visited[vertices.size()];
    for(int i = 0; i<vertices.size(); ++i){
        visited[i] = false;
    }
    
    stack <int> traversalStack;
    int sourceInt = findVertexIndex(source);
    traversalStack.push(sourceInt);
    visited[sourceInt] = true;
    
    while(!traversalStack.empty()){
        //Pop the front element
        sourceInt = traversalStack.top();
        traversalStack.pop();
        
	for(auto  v: vertices[sourceInt]->outgoingEdges) {
		if(v->label == end)
			return true;

		int adjVertexIndex = findVertexIndex(v->label);
		if(!visited[adjVertexIndex]){
			traversalStack.push(adjVertexIndex);
			visited[adjVertexIndex] = true;
		}
	}
    }

    return false;
}


		void write(llvm::raw_fd_ostream & o) const {

      o << "Graph:\n";
      for(auto const &v : vertices){
        o << "\t";
        Vertex * nonUnique_v = &(*v);
        (nonUnique_v)->write(o);
        o << "\n";
      }

		}

};

llvm::raw_fd_ostream & operator<<(llvm::raw_fd_ostream & o,
    const DiGraph & toPrint) {
  toPrint.write(o);
  return o;
}


