#include <iostream>
#include <string>

#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <map>

//#include <bits/stdc++.h>
#include <array>
#include <chrono>

#include<stack>
using namespace std;

enum direction { d_down, d_right, none };

#define COORD std::pair<int, int>    

//#define DEBUG

int iter = 0;

///Auxiliary functions

void display_arr(int* arr, int n) {

    cout << "arr: ";

    for (int i = 0; i < n; i++) {
        cout << arr[i] << " ";
    }

    cout << endl;

}

void print_coords(COORD start, COORD end) {

    cout << "Start:" << start.first << "," << start.second << endl;
    cout << "End:" << end.first << "," << end.second << endl;

}

int find_length(COORD start, COORD end, direction dir) {

    if (dir == d_down)
        return end.first - start.first;
    if (dir == d_right)
        return end.second - start.second;

    return -1;
}

void convert_sol(int** mat, int**& sol_mat, int m, int n) {

    sol_mat = new int* [m]; //Rows
    for (int i = 0; i < m; i++) {
        sol_mat[i] = new int[n]; //Cols
    }

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < m; j++) {
            if (mat[i][j] == -2)
                sol_mat[i][j] = -2; //Empty value cell
            else
                sol_mat[i][j] = -1; //Hint or empty cell
        }
    }
}

void print_one_matrix(int** matrix, int m, int n) {
    std::cout << "Matrix: " << std::endl;
    for (int i = 0; i < m; i++) { //rows
        for (int j = 0; j < n; j++) { //cols
            std::cout << matrix[i][j] << "\t";
        }
        std::cout << "\n";
    }
}

void sol_to_file(int** mat, int** sol_mat, int m, int n, string fname) {

    string fnamee = "visualize.kakuro";
    ofstream to_write(fnamee);

    to_write << m << " " << n << "\n";

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            if (mat[i][j] != -2)
                to_write << mat[i][j] << " ";
            else
                to_write << sol_mat[i][j] << " ";
        }
        to_write << "\n";
    }

    to_write.close();
}

void read_matrix(int**& matrix, std::ifstream& afile, int m, int n) {

    matrix = new int* [m]; //rows

    for (int i = 0; i < m; i++) {
        matrix[i] = new int[n]; //cols
    }

    int val;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            afile >> val;
            matrix[i][j] = val;
        }
    }
}

///Auxiliary functions

struct sum {
    COORD start;
    COORD end;

    int hint;
    int dir;
    int length;
    int* arr;

    void print_sum() {
        cout << "############################" << endl;
        cout << "Creating sum with: " << endl;
        print_coords(start, end);
        cout << "Hint: " << hint << endl;
        cout << "Direction: " << dir << endl;
        cout << "Length: " << length << endl;
        cout << "############################" << endl;
    }

    sum(COORD _start, COORD _end, int _hint, direction _dir) :
        start(_start), end(_end), hint(_hint), dir(_dir)
    {
        length = find_length(_start, _end, _dir);
        arr = new int[length];
#ifdef DEBUG
        cout << "############################" << endl;
        cout << "Creating sum with: " << endl;
        print_coords(start, end);
        cout << "Hint: " << hint << endl;
        cout << "Direction: " << dir << endl;
        cout << "Length: " << length << endl;
        cout << "############################" << endl;
#endif
    }

    //~sum(){
    //delete arr;
    //}
};


COORD find_end(int** matrix, int m, int n, int i, int j, direction dir) { //0 down 1 right

    if (dir == d_right) {
        for (int jj = j + 1; jj < n; jj++) {
            if (matrix[i][jj] != -2 || jj == n - 1) {
                if (matrix[i][jj] == -2 && jj == n - 1)
                    jj++;
                COORD END = COORD(i, jj);
                return END;
            }
        }
    }

    if (dir == d_down) {
        for (int ii = i + 1; ii < m; ii++) {
            if (matrix[ii][j] != -2 || ii == m - 1) {
                if (matrix[ii][j] == -2 && ii == m - 1)
                    ii++;
                COORD END = COORD(ii, j);
                return END;
            }
        }
    }

}


vector<sum> get_sums(int** matrix, int m, int n) {

    vector<sum> sums;

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            int val = matrix[i][j];
            if (val != -1 && val != -2) {
                int hint = val;
                hint = hint / 10;

                if ((hint % 100) == 0) {
                    hint = (int)(hint / 100);
                    COORD START = COORD(i, j + 1);
                    COORD END = find_end(matrix, m, n, i, j, d_right);
                    sum _sum = sum(START, END, hint, d_right);
                    sums.push_back(_sum);
                }

                else {
                    int div = (int)(hint / 100);
                    int rem = (int)(hint % 100);

                    if (div == 0 && rem != 0) {
                        COORD START = COORD(i + 1, j);
                        COORD END = find_end(matrix, m, n, i, j, d_down);
                        sum _sum = sum(START, END, rem, d_down);
                        sums.push_back(_sum);
                    }

                    if (div != 0 && rem != 0) {
                        COORD START1 = COORD(i + 1, j);
                        COORD START2 = COORD(i, j + 1);
                        COORD END1 = find_end(matrix, m, n, i, j, d_down);
                        COORD END2 = find_end(matrix, m, n, i, j, d_right);
                        sum _sum1 = sum(START1, END1, rem, d_down);
                        sum _sum2 = sum(START2, END2, div, d_right);
                        sums.push_back(_sum1);
                        sums.push_back(_sum2);
                    }
                }
            }


        }
    }
    return sums;
}

bool is_valid(int** sol_mat, int RowI, int ColI, int target, bool is_row, int m, int n, int digit) {
    bool Unfilled_cells = false;
    // check if the sum of the numbers in the row or column equals the target value
    int sumt = 0;
    int IterNum = 0;
    int i = 1;

    while (i < (m)) {
        IterNum = sol_mat[is_row ? RowI : i][is_row ? i : ColI];
        if (IterNum == -2 || IterNum == -1) {

            if (IterNum == -2) {
                Unfilled_cells = true;
            }
            IterNum = 0;
        }
        sumt += IterNum;
        i++;
    }
    if (Unfilled_cells) {
        if (sumt >= target || sumt == 0) {
            //cout << "unfilleds cells and above target" << sumt << "/" << target << endl;
            return false;
        }
    }
    else {
        if (sumt != target || sumt == 0) {
            //cout << "filleds cells and below target "<<sumt<<"/"<<target<<RowI<<"ggg" << ColI << is_row <<" digit is"<<digit << endl;
            return false;
        }
    }

    //print_one_matrix(sol_mat, m, n);
    //cout << endl;
    //cout <<sumt<<"/" << target;
    // check if there are any repeated digits in the row or column
    // check if there are any repeated digits in the row or column
    bool has_repeated_digits = false;

    if (is_row) {
        for (int j = 0; j < n; j++) {
            if (sol_mat[RowI][j] > 0 && j != ColI) {
                int digitt = sol_mat[RowI][j];
                if (digitt == digit) {
                    has_repeated_digits = true;
                    break;
                }
            }
        }
    }
    else {
        for (int i = 0; i < m; i++) {
            if (sol_mat[i][ColI] > 0 && i != RowI) {
                int digitt = sol_mat[i][ColI];
                if (digitt == digit) {
                    has_repeated_digits = true;
                    break;
                }
            }
        }
    }

    if (has_repeated_digits) {
        return false;
    }


    // if both checks pass, then the row or column is valid
    //cout << "dogru mu"<< is_row;
    return true;
}


bool solution(int** mat, int** sol_mat, vector<int> sums, int m, int n, int RowI, int ColI,int lenEmptySpaces) {
    //pStack<pair<rowi,coli>digit>
    //while stack is not empty continue adding elements to the back and pop() if not meet the criteria. if execution ends and stack is empty no solution
    // if execution ends and the number of elements in the stack is equal to number of empty spaces in the kakuro it is a solution
    // calculate the next row and column indices
    int next_row = ColI == m - 1 ? RowI + 1 : RowI;
    int next_col = (ColI + 1) % n;
    stack<pair<pair<int, int>, int>> pStack;
    pair<int, int> second;
    pair<pair<int, int>, int> first;
    second.first = RowI; second.second = ColI;
    first.first = second; first.second = 1;
    pStack.push(first);
    while(!pStack.empty()||pStack.size()!=lenEmptySpaces)
    {
        
        if(is_valid(sol_mat, pStack.top().first.first, pStack.top().first.second, sums[sums.size() / 2 + RowI - 1], true,  m, n, pStack.top().second)&&is_valid(sol_mat, pStack.top().first.first, pStack.top().first.second, sums[ColI - 1], false,  m, n, pStack.top().second))
        {
            pair<int, int> second;
            pair<pair<int, int>, int> first;
            second.first = RowI; second.second = ColI;
            first.first = second; first.second = 1;
            pStack.push(first);
            sol_mat[RowI][ColI] = pStack.top().second;
            int RowI = ColI == m - 1 ? RowI + 1 : RowI;
            int ColI = (ColI + 1) % n;
        }
        else
        {
            if(pStack.top().second!=9)
            {
                pStack.top().second++;
                sol_mat[RowI][ColI] = pStack.top().second;
            }
            else
            {
                pStack.pop();
                int RowI =pStack.top().first.first;
                int ColI = pStack.top().first.second;
            }
        }   
    }
    if (pStack.empty()) {
        return false;
    }
    return true;

    if (RowI == m) {
        return true;
    }

   
    
    

    // if the current cell is already filled, move on to the next one
    //if (sol_mat[RowI][ColI] != -2) {
    //    return solution(mat, sol_mat, sums, m, n, next_row, next_col);
    //}

    //// try filling the current cell with each possible digit from 1 to 9
    //for (int digit = 1; digit <= 9; digit++) {
    //    // check if the digit is valid in the current row and column
    //    sol_mat[RowI][ColI] = digit;
    //    bool is_row_valid = is_valid(sol_mat, RowI, ColI, sums[sums.size() / 2 + RowI - 1], true, m, n, digit);
    //    bool is_col_valid = is_valid(sol_mat, RowI, ColI, sums[ColI - 1], false, m, n, digit);
    //    if (is_row_valid && is_col_valid) {
    //        // if the digit is valid, fill the current cell and move on to the next one
    //        if (solution(mat, sol_mat, sums, m, n, next_row, next_col,lenEmptySpaces)) {
    //            return true;
    //        }
    //        // if the current digit doesn't lead to a solution, backtrack and try the next one
    //        sol_mat[RowI][ColI] = -2;
    //    }
    //    else {
    //        sol_mat[RowI][ColI] = -2;
    //    }
    //}

    //// if no digit leads to a solution, backtrack to the previous cell
    //return false;
}

//TO DO: Write the solution
//You can use any algorithm and data type
//Write your solution to file in main function using        () after solving it
//int counter = 0, Number_tobeadded;
//vector<int> Numsvec;

int main(int argc, char** argv) {

    //std::string filename(argv[1]);
    string filename = "C:\\Users\\ytufek\\Desktop\\personal\\okul\\4.2\\cs406\\CS406_531_HW2\\board4_1.kakuro";
    std::ifstream file;
    file.open(filename.c_str());
    map<int, int> Hash;
    int m, n;
    file >> m;
    file >> n;

    int** mat;
    read_matrix(mat, file, m, n);
    print_one_matrix(mat, m, n);

    int** sol_mat;

    convert_sol(mat, sol_mat, m, n);
    //print_one_matrix(sol_mat, m, n);

    vector<sum> sums = get_sums(mat, m, n);
    vector<int> sorted_mat;
    vector<int> Hints;
    for (int j = 0; j < sums.size(); j++) {
        //cout << sums[j].hint << sums[j].start.first << sums[j].start.second << endl;
        if (j < sums.size() / 2) {
            sorted_mat.push_back(sums[j].start.second);//indexes of columns
            Hash.insert(make_pair(sums[j].start.second, sums[j].hint));//map the values with col index and hint
        }

    }
    sort(sorted_mat.begin(), sorted_mat.end());
    for (int r = 0; r < sums.size(); r++) {
        if (r < sums.size() / 2) {
            Hints.push_back(Hash[sorted_mat[r]]);
        }
        else {
            Hints.push_back(sums[r].hint);
        }

    }

    auto start = std::chrono::high_resolution_clock::now();
    solution(mat, sol_mat, Hints, m, n, 0, 0,sums.size());
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << elapsed.count();
    print_one_matrix(sol_mat, m, n);
    sol_to_file(mat, sol_mat, m, n, "solution.kakuro");

    for (int i = 0; i < n; i++) {
        delete mat[i];
        delete sol_mat[i];
    }

    delete mat;
    delete sol_mat;

    return 0;
}
