#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <set>
#include <map>
#include <cmath>
#include <algorithm>

using namespace std;

struct Implicant {
    string mask; // e.g., "0-101"
    set<int> minterms;
    bool used = false;

    bool operator<(const Implicant& other) const {
        if (mask != other.mask) return mask < other.mask;
        return minterms < other.minterms;
    }
};

// Return standard Gray Code sequence based on bit length (1 to 3 bits)
vector<string> getGrayCode(int bits) {
    if (bits == 1) return {"0", "1"};
    if (bits == 2) return {"00", "01", "11", "10"};
    if (bits == 3) return {"000", "001", "011", "010", "110", "111", "101", "100"};
    return {};
}

// Convert mask representation (e.g., "0-10") to readable Boolean expression
string maskToExpression(const string& mask, const vector<char>& varNames) {
    string expr = "";
    for (size_t i = 0; i < mask.size(); ++i) {
        if (mask[i] == '1') {
            expr += varNames[i];
        } else if (mask[i] == '0') {
            expr += varNames[i];
            expr += '\'';
        }
    }
    return expr.empty() ? "1" : expr;
}

bool canCombine(const Implicant& a, const Implicant& b, string& combinedMask) {
    int diffCount = 0;
    combinedMask = a.mask;
    for (size_t i = 0; i < a.mask.size(); ++i) {
        if (a.mask[i] != b.mask[i]) {
            diffCount++;
            combinedMask[i] = '-';
        }
    }
    return diffCount == 1;
}

// Quine-McCluskey Algorithm
vector<Implicant> getPrimeImplicants(const vector<int>& minterms, int numVars) {
    vector<Implicant> currentGroup;
    for (int m : minterms) {
        string bin = "";
        for (int i = numVars - 1; i >= 0; --i) {
            bin += ((m >> i) & 1) ? '1' : '0';
        }
        currentGroup.push_back({bin, {m}, false});
    }

    set<Implicant> primeImplicants;

    while (!currentGroup.empty()) {
        vector<Implicant> nextGroup;
        set<string> addedMasks;

        for (size_t i = 0; i < currentGroup.size(); ++i) {
            for (size_t j = i + 1; j < currentGroup.size(); ++j) {
                string combinedMask;
                if (canCombine(currentGroup[i], currentGroup[j], combinedMask)) {
                    currentGroup[i].used = true;
                    currentGroup[j].used = true;

                    if (addedMasks.find(combinedMask) == addedMasks.end()) {
                        addedMasks.insert(combinedMask);
                        set<int> combinedMinterms = currentGroup[i].minterms;
                        combinedMinterms.insert(currentGroup[j].minterms.begin(), currentGroup[j].minterms.end());
                        nextGroup.push_back({combinedMask, combinedMinterms, false});
                    }
                }
            }
        }

        for (const auto& imp : currentGroup) {
            if (!imp.used) {
                primeImplicants.insert(imp);
            }
        }
        currentGroup = nextGroup;
    }

    return vector<Implicant>(primeImplicants.begin(), primeImplicants.end());
}

// Petrick's Method to find all minimal covering combinations
vector<vector<Implicant>> solvePetrick(const vector<Implicant>& PIs, const vector<int>& minterms) {
    map<int, vector<int>> coverage;
    for (size_t i = 0; i < PIs.size(); ++i) {
        for (int m : PIs[i].minterms) {
            coverage[m].push_back(i);
        }
    }

    set<set<int>> pos;
    pos.insert({});

    for (int m : minterms) {
        set<set<int>> newPos;
        for (const auto& term : pos) {
            for (int piIdx : coverage[m]) {
                set<int> newTerm = term;
                newTerm.insert(piIdx);
                newPos.insert(newTerm);
            }
        }
        pos = newPos;
    }

    size_t minSize = 1e9;
    for (const auto& cover : pos) {
        if (cover.size() < minSize) {
            minSize = cover.size();
        }
    }

    vector<vector<Implicant>> minCovers;
    for (const auto& cover : pos) {
        if (cover.size() == minSize) {
            vector<Implicant> solution;
            for (int idx : cover) {
                solution.push_back(PIs[idx]);
            }
            minCovers.push_back(solution);
        }
    }

    return minCovers;
}

int main() {
    string filename = "kmap.txt";
    ifstream infile(filename);

    if (!infile.is_open()) {
        cerr << "Error: Could not open " << filename << endl;
        return 1;
    }

    vector<vector<int>> matrix;
    string line;

    // Read the matrix line by line to determine order dynamically
    while (getline(infile, line)) {
        stringstream ss(line);
        vector<int> row;
        int val;
        while (ss >> val) {
            row.push_back(val);
        }
        if (!row.empty()) {
            matrix.push_back(row);
        }
    }
    infile.close();

    int R = matrix.size();
    if (R == 0) {
        cerr << "Error: Empty input file." << endl;
        return 1;
    }
    int C = matrix[0].size();

    // Variable count determination
    int rowBits = round(log2(R));
    int colBits = round(log2(C));
    int numVars = rowBits + colBits;

    if (numVars < 1 || numVars > 5) {
        cerr << "Error: Matrix dimensions (" << R << "x" << C 
             << ") correspond to " << numVars << " variables. Only 1 to 5 variables supported." << endl;
        return 1;
    }

    vector<char> allVars = {'a', 'b', 'c', 'd', 'e'};
    vector<char> activeVars(allVars.begin(), allVars.begin() + numVars);

    vector<string> rowGray = getGrayCode(rowBits);
    vector<string> colGray = getGrayCode(colBits);

    // Extract minterms mapping through Gray code
    vector<int> minterms;
    for (int r = 0; r < R; ++r) {
        for (int c = 0; c < C; ++c) {
            if (matrix[r][c] == 1) {
                string bin = (rowBits > 0 ? rowGray[r] : "") + (colBits > 0 ? colGray[c] : "");
                int mintermIndex = stoi(bin, nullptr, 2);
                minterms.push_back(mintermIndex);
            }
        }
    }

    cout << "Detected K-Map: " << R << "x" << C << " (" << numVars << " variables: ";
    for (size_t i = 0; i < activeVars.size(); ++i) {
        cout << activeVars[i] << (i + 1 < activeVars.size() ? ", " : "");
    }
    cout << ")\n" << endl;

    int totalMinterms = 1 << numVars;
    if (minterms.empty()) {
        cout << "Minimized Expression: 0" << endl;
        return 0;
    }
    if ((int)minterms.size() == totalMinterms) {
        cout << "Minimized Expression: 1" << endl;
        return 0;
    }

    // Solve for minimal expressions
    vector<Implicant> PIs = getPrimeImplicants(minterms, numVars);
    vector<vector<Implicant>> minimalCovers = solvePetrick(PIs, minterms);

    cout << "Possible Minimized Boolean Expression(s):" << endl;
    for (const auto& cover : minimalCovers) {
        string result = "";
        for (size_t i = 0; i < cover.size(); ++i) {
            result += maskToExpression(cover[i].mask, activeVars);
            if (i + 1 < cover.size()) {
                result += " + ";
            }
        }
        cout << result << endl;
    }

    return 0;
}
