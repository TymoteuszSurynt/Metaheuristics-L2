#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <random>
#include <fstream>
#include <zconf.h>

using namespace std;
typedef mt19937 rng_type;

rng_type rng;
//struct
struct city {
    double x;//coordinates
    double y;//coordinates
    double *distance; //table of distances to any other city
    bool visited; //for the first path
};

//Utils
//print path
void printGraph(int *result, vector<struct city> cities, int size) {
    int check=-1;
    for (int i = 0; i < size; i++) {
        if(check!=-1) {
            cerr << result[i] + 1 << endl;
        }else if(check ==-1){
            if(result[i]==0){
                check=i;
                cerr << result[i] + 1 << endl;
            }
        }
    }
    for (int j = 0; j < check; ++j) {
        cerr << result[j] + 1 << endl;
    }
    cerr<<1<<endl;
}
//sum of distances in the result
long double sumOfDistances(int *result, vector<struct city> &cities, int size) {
    long double sum = 0;
    for (int i = 0; i < size - 1; i++) {
        sum += cities.at((unsigned long long int) result[i]).distance[result[i + 1]];
    }
    sum+=cities.at((unsigned long long int) result[size - 1]).distance[result[0]];
    return sum;
}

//Counting distance
double distance(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x1 - x2, 2.0) + pow(y1 - y2, 2.0));
}

//function for first path
int greedyDistance(struct city city, int size, vector<struct city> &cities) {
    int min = -1;
    unsigned i = 0;
    double minValue = -1;
    for (; i < size; i++) {
        if (city.distance[i] != 0) {
            if (!cities.at(i).visited) {
                minValue = city.distance[i];
                min = i;
                break;
            }
        }
    }
    if (minValue == -1) {
        return -1;
    }
    for (; i < size; ++i) {
        if (city.distance[i] < minValue && city.distance[i] != 0) {
            if (!cities.at(i).visited) {
                minValue = city.distance[i];
                min = i;
            }
        }
    }
    return min;
}
//function allowing to change 2 nodes in the graph
void opt2(int *result, int *result2, int size, int x, int y) {
    int i, j = 0;
    if(result[x]!=result[y]) {
        for (i = 0; i < x; i++) {
            result2[i] = result[i];
        }
        for (i = x; i < y + 1; i++) {
            result2[i] = result[y - j];
            j++;
        }
        for (i = y + 1; i < size; i++) {
            result2[i] = result[i];
        }
    }
}

int main() {
    std::ifstream in("t5.test");
    std::streambuf *cinbuf = std::cin.rdbuf(); //save old buf
    std::cin.rdbuf(in.rdbuf()); //redirect std::cin to in.txt!
    clock_t begin = clock();
    int a, cityName, h = 0, w = 1;
    double xCity, yCity;
    long double pathSum = 0,radius,sum, temp=1.0,tempmin=0.000000001,alpha=0.9999;
    string input;
    vector<struct city> cities;
    cin >> a;
    int *result = new int[a];
    for (int i = 0; i < a; i++) {
        struct city city;
        cin >> cityName;
        cin >> xCity >> yCity;
        city.x = xCity;
        city.y = yCity;
        city.visited = false;
        cities.push_back(city);
    }
    for (unsigned i = 0; i < a; i++) {
        cities.at(i).distance = new double[a];
        for (int j = 0; j < a; j++) {
            cities.at(i).distance[j] = distance(cities[i].x, cities[i].y, cities[j].x, cities[j].y);
        }
    }
    result[0]=0;
    while (true) {
        cities.at((unsigned long long int) h).visited = true;
        h = greedyDistance(cities.at((unsigned long long int) h), a, cities);
        if (h == -1) {
            pathSum += cities.at((unsigned long long int) result[w - 1]).distance[0];
            w++;
            break;
        } else {
            result[w] = h;
            pathSum += cities.at((unsigned long long int) result[w - 1]).distance[result[w]];
            w++;
        }

    }
    rng_type::result_type const seedval = (const unsigned int) time(NULL);
    rng.seed(seedval);
    int *result2 = new int[a];
    if(a>8000) {
        radius = pathSum/a;
    }else{
        radius=pathSum;
    }
    uniform_int_distribution<int> udist(0, a - 1);
    uniform_real_distribution<double> udist1(0.0, 1.0);
    while (temp>tempmin) {
        int random_number = udist(rng);
        int random_number1 = udist(rng);
        if(cities[result[random_number]].distance[result[random_number1]]<radius) {
            if (random_number != random_number1) {
                opt2(result, result2, a, result[random_number], result[random_number1]);
                sum = sumOfDistances(result2, cities, a);
                double random_number2 = udist1(rng);
                if (sum <= pathSum) {
                    for (int l = 0; l < w; ++l) {
                        result[l] = result2[l];
                    }
                    pathSum = sum;
                } else if (exp((pathSum - sum) / temp) > random_number2) {
                    for (int l = 0; l < w; ++l) {
                        result[l] = result2[l];
                    }
                    pathSum = sum;
                }
                temp = temp * alpha;
            }
        }
        if(double(clock()-begin)/ CLOCKS_PER_SEC>599){
            break;
        }
    }

    delete[] result2;
    cout<<pathSum<<endl;
    printGraph(result, cities, a);
    for (unsigned i = 0; i < a; i++) {
        delete[] cities.at(i).distance;
    }
    return 0;
}