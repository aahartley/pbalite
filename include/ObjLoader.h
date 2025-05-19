#ifndef OBJ_LOADER_H
#define OBJ_LOADER_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "Vector.h"

namespace pba
{

struct Triangle
{
  Triangle(){}
  Triangle(Vector vv1, Vector vv2, Vector vv3) : v1(vv1), v2(vv2), v3(vv3)
  {
    e1 = v2 - v1;
    e2 = v3 - v1;
    e3 = v3 - v2;
  }

  Vector v1, v2, v3, e1, e2, e3;
};
class ObjLoader {
public:
    static std::vector<Triangle> loadObj(const std::string& filename) {
        std::vector<Vector> vertices;
        std::vector<Triangle> triangles;

        std::ifstream file(filename);
        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string token;
            iss >> token;

            if (token == "v") {
                Vector vertex;
                float x,y,z;
                iss >> x >> y >> z;
                vertex.set(x,y,z);
                vertices.push_back(vertex);
            } else if (token == "f") {
                std::string v1s, v2s, v3s;
                int v1Index, v2Index, v3Index;
                iss >> v1s  >> v2s  >> v3s;
                size_t pos1 = v1s.find('/'); // Find the position of the first '/'
                size_t pos2 = v2s.find('/'); // Find the position of the first '/'
                size_t pos3 = v3s.find('/'); // Find the position of the first '/'
                if (pos1 != std::string::npos)
                { // Check if the '/' is found
                    std::string num = v1s.substr(0, pos1); // Extract from start to the position of '/'
                    v1Index = std::stoi(num);
                    num = v2s.substr(0, pos2); // Extract from start to the position of '/'
                    v2Index = std::stoi(num);        
                    num = v3s.substr(0, pos3); // Extract from start to the position of '/'
                    v3Index = std::stoi(num);
                }
                else
                {
                    v1Index = std::stoi(v1s);
                    v2Index = std::stoi(v2s);
                    v3Index = std::stoi(v3s);

                }
                Vector v1 = vertices[v1Index - 1];
                Vector v2 = vertices[v2Index - 1];
                Vector v3 = vertices[v3Index - 1];
                triangles.push_back(Triangle(v1,v2,v3));
                //std::cout << v1Index << ' ' << v2Index << ' ' << v3Index << ' ' << v1.X() << ' ' << v1.Y() << ' ' << v1.Z() << '\n';
            }
        }

        file.close();
        return triangles;
    }
    static std::vector<Vector> loadObjVert(const std::string& filename) {
        std::vector<Vector> vertices;
        std::vector<Triangle> triangles;

        std::ifstream file(filename);
        std::string line;

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string token;
            iss >> token;

            if (token == "v") {
                Vector vertex;
                float x,y,z;
                iss >> x >> y >> z;
                vertex.set(x,y,z);
                vertices.push_back(vertex);
            } else if (token == "f") {
                std::string v1s, v2s, v3s;
                int v1Index, v2Index, v3Index;
                iss >> v1s  >> v2s  >> v3s;
                size_t pos1 = v1s.find('/'); // Find the position of the first '/'
                size_t pos2 = v2s.find('/'); // Find the position of the first '/'
                size_t pos3 = v3s.find('/'); // Find the position of the first '/'
                if (pos1 != std::string::npos)
                { // Check if the '/' is found
                    std::string num = v1s.substr(0, pos1); // Extract from start to the position of '/'
                    v1Index = std::stoi(num);
                    num = v2s.substr(0, pos2); // Extract from start to the position of '/'
                    v2Index = std::stoi(num);        
                    num = v3s.substr(0, pos3); // Extract from start to the position of '/'
                    v3Index = std::stoi(num);
                }
                else
                {
                    v1Index = std::stoi(v1s);
                    v2Index = std::stoi(v2s);
                    v3Index = std::stoi(v3s);

                }
                Vector v1 = vertices[v1Index - 1];
                Vector v2 = vertices[v2Index - 1];
                Vector v3 = vertices[v3Index - 1];
                triangles.push_back(Triangle(v1,v2,v3));
                //std::cout << v1Index << ' ' << v2Index << ' ' << v3Index << ' ' << v1.X() << ' ' << v1.Y() << ' ' << v1.Z() << '\n';
            }
        }

        file.close();
        return vertices;
    }
};
    
} 
#endif // OBJ_LOADER_H