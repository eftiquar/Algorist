//
//  main.cpp
//  Uber
//
//  Created by zayan on 11/16/19.
//  Copyright © 2019 Eftiquar. All rights reserved.
//

#include <iostream>
#include <vector>
#include <map>
#include <bitset>
#include <queue>
#include <algorithm>
#include <cmath>
#include <set>
using namespace std;
using std::bitset;
using std::vector;
using std::map;
using std::string;
using std::wstring;
using std::cout;

/*
      t     e     l     l    o
   0  1     2     3     4    5
 t 1  0     1     2     3    4
 h 2  1
 e 3  2     1     2     3    4
 l 4  3     2     1     2    3
 l 5  4     3     2     1    2
 o 5  5     4     3     2    1
*/
class SquareGridWalker
{
public:
    bool IsOnLeftDiagonal(int i, int j)const
    {
        return  false;
    }
    SquareGridWalker(const vector<vector<int>>& grid):grid_(grid),
    iMax_(grid[0].size()),jMax_(iMax_)
    {
        
    }
private:
    const vector<vector<int>>& grid_;
    const size_t iMax_, jMax_;
};
unsigned char Reverse(unsigned char in)
{
    unsigned char result = 0;
    unsigned char mask = 1 << 7;
    for(int i = 0; i != 8;++i)
    {
        result >>= 1;
        result |= in & mask;
        in <<=1;
    }
    return result;
}

class LongPal
{
public:
    LongPal(string input):input_(input)
    {
        input_length_ = input_.size();
        dp_ = vector<vector<size_t>>(input_length_,vector<size_t>(input_length_));
        
    }
    size_t LongPalinDrome()
    {
        //liril
        
        for(int j =0;j!= input_length_; ++j )
        {
            for(int i = j + 1; i !=0; )
            {
                --i;
               if(input_[i] == input_[j])
               {
                    if(j - i == 1)
                    {
                            dp_[i][j] = 2;
                    }
                   else if(j == i)
                   {
                            dp_[i][j] = 1;
                   }
                   else
                   {
                       dp_[i][j] = dp_[i+1][j-1] + 2;
                   }
               }
               if(input_[i] != input_[j])
               {
                   size_t right_max = dp_[i + 1][j];
                   size_t left_max =  dp_[i][j-1];
                   auto max_new = max(right_max,left_max);
                   dp_[i][j] = max(dp_[i][j],max_new);
               }
            }
        }
        return dp_[0][input_length_-1];
    }
    string Reconstruct()
    {
        string result;
        for(size_t i = input_length_;i;)
            for(size_t j = input_length_;j;)
        {
            --i;
            --j;
            
        }
        return result;
    }
private:
    string input_;
    size_t input_length_;
    vector<vector<size_t>> dp_;
};
class KnapSack
{
public:
    KnapSack(vector<pair<int,int>> items):items_(items)
    {
        
    }
    int Compute(int capacity)
    {
        vector<int> dp(capacity + 1);
        for(int item = 0; item != items_.size();++item)
        {
            for(int kCap = capacity;kCap < items_[item].first == false; --kCap  )
            {
                dp[kCap] = max(dp[kCap],dp[kCap - items_[item].first] + items_[item].second);
            }
        }
        return dp.back();
    }
private:
    vector<pair<int,int>> items_;
};

class EditDistance
{
public:
    EditDistance(string src, string target):src_(src), target_(target)
    {
       
        srcVector_.push_back(' ');
        std::copy(src.cbegin(), src.cend(), std::back_inserter(srcVector_));
        
       
        targetVector_.push_back(' ');
        std::copy(target_.cbegin(), target_.cend(), std::back_inserter(targetVector_));
        
        std::vector<int> vectorX0;
        for(int i = 0; i !=target_.size() +1;++i)
        {
            vectorX0.push_back(i);
        }
        
        editDP_ = std::vector<std::vector<int>>(srcVector_.size(),std::vector<int>(targetVector_.size()));
        editDP_[0] =vectorX0;
        for(int i = 0; i != srcVector_.size();++i)
        {
            editDP_[i][0] = i;
        }
        
        
            for(int x = 1; x!= targetVector_.size();++x )
                for(int y = 1;y != srcVector_.size();++y)
            {
                 if(srcVector_[y] == targetVector_[x])
                     editDP_[y][x] = editDP_[y-1][x-1];
                else
                {
                    editDP_[y][x] = editDP_[y-1][x-1] + 1;
                    int delCost = editDP_[y-1][x] + 1;
                    int insertCost = editDP_[y][x-1] + 1;
                    int minCost = std::min(insertCost,delCost);
                    if(minCost < editDP_[y][x])
                        editDP_[y][x] = minCost;
                }
            }
        PrintSolution(targetVector_.size()-1,srcVector_.size()-1);
    }
private:
    void PrintSolution(size_t x, size_t y)
    {
        if(x ==0 && y ==0)
            return;
        if(x==0)
        {
            PrintSolution(x,y-1);
            cout << "delete:" << srcVector_[y] << std::endl;
        }
        else if(y ==0)
        {
            PrintSolution(x-1,y);
            cout << "insert:" << targetVector_[x] << std::endl;
        }
        else if(editDP_[y][x] == editDP_[y-1][x-1] && targetVector_[x] == srcVector_[y])
        {
            PrintSolution(x-1, y-1);
            std::cout << "match:" << targetVector_[x] << std::endl;
        }
        else if(editDP_[y][x] == editDP_[y-1][x-1] +1)
        {
            PrintSolution(x-1, y-1);
            std::cout << "sub:" << targetVector_[x] << std::endl;
        }
        else if(editDP_[y][x] == editDP_[y-1][x] + 1)
        {
            PrintSolution(x, y-1);
            cout << "delete:" << srcVector_[y] << std::endl;
        }
        else if(editDP_[y][x] == editDP_[y][x-1] + 1)
        {
            PrintSolution(x-1, y);
            cout << "insert:" << targetVector_[x] << std::endl;
        }
    }
    string src_, target_;
    std::vector<std::vector<int>> editDP_;
    std::vector<char> srcVector_,targetVector_;
};

class Sudoku
{
public:
    Sudoku(vector<vector<size_t>>& sudoku_grid):sudoku_grid_(sudoku_grid)
    {
        if(Solve(begin,begin))
            Print();
    }
    void Print()
    {
        for(int x = begin; x!= end; ++x)
        {
            std::ostream_iterator<int> console(std::cout,",");
            std::copy(sudoku_grid_[x].cbegin(), sudoku_grid_[x].cend(), console);
            std::cout << "\n";
        }
    }
    bool Solve(int x,int y)
    {
        if(x == end )
        {
            if(++y == end)
                return true;
            x = begin;
        }
        if(sudoku_grid_[x][y] != kBlankValue)
            return Solve(x+1,y);
        
        for(int val = begin+1; val != end+1;++val)
        {
            if(CanPlaceVal(x, y, val))
            {
                sudoku_grid_[x][y] = val;
                if(Solve(x+1,y))
                    return true;
                sudoku_grid_[x][y] = kBlankValue;
                
            }
        }
       
        return false;
    }
private:
    bool CanPlaceVal(int x, int y, int val)
    {
        return IsXConflict(x, val) || IsYConflict(y, val) || IsCellConflict(x, y, val) ? false : true;
    }
    bool IsXConflict(int x,  int val)const
    {
        return std::find(sudoku_grid_[x].begin(),sudoku_grid_[x].end(),val) !=sudoku_grid_[x].end();
    }
    bool IsYConflict( int y, int val)const
    {
        for(int x = 0;x != end; ++x)
        {
            if(sudoku_grid_[x][y] == val)
                return true;
        }
        return false;
    }
    bool IsCellConflict(int x, int y, int val)const
    {
        int cellxbegin = (x/3) * 3;
        int cellybegin = (y/3) * 3;
        
        for(int xCursor = cellxbegin; xCursor != cellxbegin + 3;++xCursor)
            for(int yCursor = cellybegin;yCursor != cellybegin +3; ++yCursor)
            {
                if(sudoku_grid_[xCursor][yCursor] == val)
                    return true;
            }
        return false;
    }
    const int begin = 0;
    const int end = 9;
    const int kBlankValue = 0;
    vector<vector<size_t>>& sudoku_grid_;
    
    };
int alternatingCharacters(string s) {

    int result = 0;
    auto first = s.begin();
    for( auto next = first+1; next != s.end();++next,++first)
    {
        if(*next == *first)
        {
            ++result;
        }
        
    }
    return result;
}
void revert(vector<int>::iterator first,vector<int>::iterator last )
{
    while(first != last && first != --last)
    {
        iter_swap(first++,last);
    }
}

//123 >> 231 >> 312 >> 123
vector<int> rotLeft(vector<int> a, int d)
{
    d %= a.size();
    if(d)
    {
        revert(a.begin(),a.begin()+d);
        revert(a.begin()+d,a.end());
        revert(a.begin(),a.end());
    }
    return a;
}
int makeAnagram(string a, string b)
{
    map<char,int> a_map;
    map<char,int> b_map;
    for(const auto ch: a)
    {
        ++a_map[ch];
    }
    for(const auto ch: b)
    {
        ++b_map[ch];
    }
    for( auto& a_entry: a_map)
    {
        auto b_entry =b_map.find(a_entry.first);
        if(b_entry != b_map.end())
        {
            if(a_entry.second < b_entry->second)
            {
                b_entry->second -= a_entry.second;
                a_entry.second = 0;
            }
            else if(a_entry.second == b_entry->second)
            {
                a_entry.second = 0;
                b_entry->second = 0;
            }
            else
            {
                a_entry.second -= b_entry->second;
                b_entry->second =0;
            }
        }
    }
    int result = 0;
    for(const auto& a_entry : a_map)
    {
        result += a_entry.second;
    }
    for(const auto& b_entry: b_map)
    {
        result += b_entry.second;
    }
    return result;
}
int maximumToys(vector<int> prices, int k)
{
    sort(prices.begin(),prices.end());
    int consumed = 0;
    int i = 0;
    for(; i != prices.size();++i)
    {
        if(consumed + prices[i] > k == false)
        {
            consumed += prices[i];
        }
        else
            break;
    }
    return i;
}

class JumpGame
{
public:
    JumpGame(vector<int> game):game_(game),
    end_game_(game.size())
    {
        
    }
    vector<size_t> Solve()
    {
        return solution_;
    }
    int BFS()
    {
        //[2,3,1,1,4]
        
        queue<int> next_steps;
        next_steps.push(0);
        vector<int> dp(end_game_,numeric_limits<int>::max());
        dp[0] = 0;
        while(!next_steps.empty())
        {
            auto this_step = next_steps.front();
            next_steps.pop();
            int max_reach = this_step + game_[this_step];
            if(max_reach < end_game_ -1  == false)
            {
                return dp[this_step] + 1;
            }
            auto end_step = min(end_game_,max_reach + 1);
            for(int i = this_step + 1; i < end_step; ++i)
            {
                dp[i] = dp[this_step] + 1;
                next_steps.push(i);
            }
        }
        return -1;
    }
    int MinJumps()
    {
    
        
        //[2,3,1,1,4]
        /*
         step 0 - max 2, step count 1
         step 1 - max 4, step count 2
         step 2 - max 3, step count 2
         
         */
        vector<int> dp(end_game_,numeric_limits<int>::max());
        
        dp[0] = 0;
        int max_jump = std::numeric_limits<int>::min();
        int min_steps = std::numeric_limits<int>::max();
        for(int i = 0; i != end_game_; ++ i)
        {
            max_jump = max(max_jump, i + game_[i]);
            min_steps = min(min_steps,dp[i] + 1);
        }
        /*
         dp[0] = 0;
        for(int i = 1; i != end_game_; ++ i)
        {
            for(int j = i; j;)
            {
                --j;
                if(j + game_[j] < i == false)
                {
                    dp[i] = min(dp[i],dp[j] + 1);
                }
                //
            }
        }*/
        //for(int )
        return dp.back();
    }
private:
    bool BackTrack(size_t cursor)
    {
        if(cursor < end_game_ == false)
        {
            return true;
        }
        size_t max_jump = game_[cursor] + cursor;
        for(size_t try_jump = max_jump;try_jump != cursor; --try_jump)
        {
           if(BackTrack(try_jump))
           {
               solution_.push_back(try_jump);
               return true;
           }
        }
        return false;
    }
    vector<size_t> solution_;
    
    vector<int> game_;
    int end_game_;
};
class ScoreWays
{
public:
    ScoreWays(vector<int> scores, int target):target_(target), scores_(move(scores))
    {
        
    }
    int Combinations()
    {
        vector<int> dp(target_ + 1);
        dp[0] = 1;
        for(auto score : scores_)
        {
            for(int i = 1; i != target_; ++i)
            {
                if(score > i == false)
                {
                    dp[i] += dp[i-score];
                }
                //dp[i] = score <= i ? dp[i-score]
            }
        }
        return dp.back();
    }
    int Permutations()
    {
        vector<int> dp(target_ + 1);
        dp[0] = 1;
        for(int i = 1; i != target_; ++i)
        {
            for(auto score : scores_)
            {
                if(score > i == false)
                {
                    dp[i] += dp[i-score];
                }
                //dp[i] = score <= i ? dp[i-score]
            }
        }
        return dp.back();
        
    }
private:
    int target_;
    vector<int> scores_;
};
//1,2,4,2,3,6,1,5
class EqualSumBackTrack
{
public:
    EqualSumBackTrack(vector<int>input, int target):input_(input),target_(target),used_(input.size())
    {
        CanProceed();//assume input is valid
    }
    bool DoBackTrack(int index_cursor)
    {
      if(index_cursor == input_end_)
          return running_target_ == 0;
        
        if(used_[index_cursor])
         return   DoBackTrack(++index_cursor);
        if(running_target_ + input_[index_cursor] > target_ == false)
        {
            running_target_ += input_[index_cursor];
            used_[index_cursor] = true;
            running_solution_.push_back(input_[index_cursor]);
            if(running_target_ == target_)
            {
                running_target_ = 0;
                if(DoBackTrack(0))
                {
                    ++running_partition_;
                    return true;
                }
                running_target_ = target_ - input_[index_cursor];
                used_[index_cursor] = false;
                running_solution_.pop_back();
                return DoBackTrack(index_cursor + 1);
            }
            if(DoBackTrack(index_cursor + 1))
            {
                return true;
            }
            else
            {
                running_target_ -= input_[index_cursor];
                used_[index_cursor] = false;
                running_solution_.pop_back();
                return DoBackTrack(index_cursor + 1);
            }
        }
        return running_partition_ == partitions_;
    }
private:
    bool CanProceed()
    {
        int sum = 0;
        for(auto elem : input_)
        {
            sum += elem;
        }
        if(sum % target_ != 0)
        return false;
        partitions_ = sum / target_;
        input_end_ = input_.size();
        
        return true;
    }
    size_t input_end_ = 0;
    int partitions_ = 0;
    int target_;
    vector<int> input_;
    vector<bool> used_;
    int running_target_=0;
    int running_partition_=0;
    vector<vector<int>> solution_;
    vector<int> running_solution_;
};
long substrCount(int n, string s)
{

    return 0;
}

/*
 a b c
 d e
 f g
 
 
 a -> a  b  c
      4  1  0
 b -> b  c  d
      4  1  0
 c -> c  d  e
      4  1  0
 d -> d  e  f
      4  1  0
 e -> e  f  g
      4  1  0
 f -> f  g
      4  1
 g -> g
      4
 */
/*
   target sum
   0 1 2 3 4 5 6 7 8 9 10
 0 1 0 0 0 0 0 0 0 0 0
 1 1 1 2 3 4 5 6 7 8 9
 2 1 1 1 2 2 3 3 4 4
 3 0 1 1 2
 4 0
 5 0
 6
 7
 8
 9
 10
 
 */
class MinCoins
{
public:
    MinCoins(const vector<int>& coins):coins_(coins)
    {
        
    }
    
    int DoMinFast(int sum)
    {
        vector<int> dp(sum+1);
        for(int i= 1; i != sum +1; i ++)
        {
            dp[i] = i;
        }
        for(int sumCursor = 2; sumCursor != sum +1; ++sumCursor)
        {
            for(int coin = 1; coin != coins_.size(); ++coin )
            {
                int cost_w_coin = sumCursor< coins_[coin] ? dp[sumCursor] : dp[sumCursor - coins_[coin]] + 1;
                cost_w_coin < dp[sumCursor] ? dp[sumCursor] = cost_w_coin :0;
            }
        }
        return dp.back();
    }
    
    int DoMinCoin(int sum)
    {
        
        vector<vector<int>> mincoins(coins_.size(),vector<int>(sum +1));
        for(int i = 0; i != sum +1;++i)
        {
            mincoins[0][i] = i;
        }
        for(int coin = 1;coin != coins_.size();++coin)
        {
            mincoins[coin][0] = 0;
            for(int thissum =1; thissum != sum +1;++thissum)
            {
                if( thissum < coins_[coin])
                {
                    mincoins[coin][thissum] = mincoins[coin-1][thissum];
                }
                else
                {
                    int with_coin = mincoins[coin][thissum-coins_[coin]] + 1;
                    int without_coin = mincoins[coin-1][thissum];
                    mincoins[coin][thissum] = with_coin < without_coin ? with_coin : without_coin;
                }
            }
        }
        mincoins_ = std::move(mincoins);
        return mincoins_.back()[sum];
    }
    vector<int> GetMinCoins(int sum)
    {
        vector<int> result;
        int i = coins_.size()-1;
        int j = sum;
        while(j)
        {
            if((i == 0 )|| mincoins_[i][j] != mincoins_[i-1][j])
            {
                j -= coins_[i];
                result.push_back(coins_[i]);
            }
            else
            {
                --i;
            }
        }
        return  result;
        
    }
    int DoWays(int sum)
    {
        vector<vector<int>> mincoins(coins_.size(),vector<int>(sum +1));
        //mincoins[0][0] = 1;
        for(int i = 0; i != sum +1;++i)
        {
            mincoins[0][i] = 1;
        }
        for(int coin = 1; coin != coins_.size();++coin)
        {
            mincoins[coin][0] = 1;
            for(int targ = 1; targ != sum + 1; ++targ)
            {
                int ways_wo_coin = mincoins[coin-1][targ] ;
                if(targ < coins_[coin])
                {
                    mincoins[coin][targ] = ways_wo_coin ;
                }
                else
                {
                    int ways_with_coin = mincoins[coin-1][targ - coins_[coin]];
                    mincoins[coin][targ] = ways_wo_coin + ways_with_coin;
                }
            }
        }
        mincoins_ = std::move(mincoins);
        return mincoins_.back()[sum];
    }

    const vector<int> coins_;
    vector<vector<int>> mincoins_;
    
};
class MaxSquare
{
  
    
};
class SuffixArray
{
public:
    SuffixArray(const wchar_t* input_string)
    {
        wstring input{input_string};
        auto len = input.size();
        
        for(int i =0; i != len; ++i)
        {
            suffixes_.push_back(input_string+i);
        }
    
        std::sort(suffixes_.begin(), suffixes_.end(),
                  [](auto first,auto second)
        {
            return !(0 < wcscmp(first,second));
        });
        
    }
private:
    vector<const wchar_t*> suffixes_;
    
};
using std::queue;
using std::pair;
class Regional
{
public:
  Regional(const vector<vector<int>>& region):region_(region)
    {
        i_max = region_.size();
        j_max = region[0].size();
        visted_ = vector<bool>(i_max*j_max);
    }
    int GetMax()
    {
        auto pos_max = visted_.size();
        int max = 0;
        for(size_t pos =0; pos != pos_max; ++ pos )
        {
            if(int this_max = DoBFS(pos))
            {
                if(max < this_max)
                    max = this_max;
            }
                
        }
        return max;
    }
private:
    int DoBFS(int position)
    {
        auto [i,j] = GetIJ(position);
        if(IsVisited(i,j))
            return 0;
        
        int result = 0;
        
        queue<pair<int,int>> bfs_queue;
        bfs_queue.push({i,j});
        SetVisited(i, j);
        while(!bfs_queue.empty())
        {
            auto popped = bfs_queue.front();
            bfs_queue.pop();
            ++result;
            
            for(auto [i_next,j_next] : GetNeighbors(popped.first,popped.second))
            {
                if(!IsVisited(i_next, j_next))
                {
                    bfs_queue.push({i_next,j_next});
                    SetVisited(i_next, j_next);
                }
            }
            
        }
        return  region_[i][j] ?  result:0;
    }
    bool IsVisited(int i, int j)
    {
        return visted_[i*j_max +j];
    }
    void SetVisited(int i, int j)
    {
        visted_[i*j_max +j] = true;
    }
    vector<std::pair<int,int>> GetNeighbors(int i, int j)
    {
        vector<std::pair<int,int>> result;
        for(const auto& [i_offset, j_offset] : neighbour_offsets_)
        {
            int new_i = i + i_offset;
            int new_j = j + j_offset;
            
            if(new_i < 0 || new_i > i_max -1)
              continue;
            if(new_j < 0 || new_j > j_max -1)
                continue;
            if(region_[i][j] == region_[new_i][new_j])
            {
                result.push_back({new_i,new_j});
            }
        }
        return  result;
        
    }
    vector<std::pair<int,int>> neighbour_offsets_ =
    {
        {-1,-1},{-1,0},{-1,1},
            {0,-1},{0,1},
        {1,-1},{1,0},{1,1}
    };
    std::pair<int,int> GetIJ(int position)
    {
        return {position/(j_max) ,position%(j_max)};
    }
    vector<bool> visted_;
    size_t i_max = 0;
    size_t j_max = 0;
    const vector<vector<int>>& region_;
};
int stepPerms(int n)
{
if (n <3)
    return n;
    int step_n_3 = 1;
    int step_n_2 = 1;
    int step_n_1 = 2;
    int result = 0;
    for (int step = 3; step != n + 1; ++step)
    {
        result = step_n_3 + step_n_2 + step_n_1;
        step_n_3 = step_n_2;
        step_n_2 = step_n_1;
        step_n_1 = result;
    }
    return result;
}
/*
 case 1
 myfreq > last
 if(diff > 1)
 1,3
 3,1
 
 my freq < last
 
 */
string isValid(string s) {

    map<char,int> char_freq;
    for(auto ch: s)
    {
        ++char_freq[ch];
    }
    int unique_chars = char_freq.size();
    int mod = s.size()% unique_chars;
    if(mod <2)
        return "YES";
    unique_chars -=1;
    if((s.size() -1) % unique_chars ==0)
        return "YES";
    int ones = 0;
    for(auto entry: char_freq)
    {
        if(entry.second ==1)
            ++ones;
    }
    
    return ones == 1 ? "YES" : "NO";
}
int lcs (string lhs, string rhs)
{
    vector<int> lcs_dp(rhs.size() + 1);
    
    for(int i = 0; i != lhs.size(); ++i)
    {
        int next_left_top = 0;
        for(int j = 1; j != lcs_dp.size(); ++j)
        {
            int leftTop = next_left_top;
            next_left_top = lcs_dp[j];
            if(lhs[i] == rhs[j-1])
            {
                lcs_dp[j] = leftTop +1;
            }
            else
            {
                lcs_dp[j] = std::max(lcs_dp[j-1],lcs_dp[j]);
            }
        }
    }
    return lcs_dp.back();
}

int poisonousPlants(vector<int> p)
{
    if(p.size() == 0)
        return 0;
    //3,6,2,7,5
    //for(int i)
    int result = 0;
    int minimum = std::numeric_limits<int>::max();
    for(auto elem: p)
    {
        if(elem < minimum)
            minimum = elem;
        else if(elem > minimum)
        {
            ++result;
        }
    }
    
    
    return result;
}
/*
 1,1,1,2,1,1,1,1
 
 */
using std::move;
class EnumSubstrings
{
    using Iter = string::iterator;
public:
    EnumSubstrings(string input):input_(move(input))
    {
        
    }
    void Enum()
    {
        EnumHelper(input_.begin(),input_.end());
    }

    int CountAnagrams()
    {
        map<string,vector<string>> anagrams;
        for(const auto& substr:substrings_)
        {
          if(substr.size() == 1)
              anagrams[substr].push_back(substr);
          else
          {
              auto key = substr;
              std::sort(key.begin(), key.end());
              anagrams[key].push_back(substr);
          }
        }
        int result = 0;
        for(auto [key,entries]:anagrams)
        {
            if(entries.size() >1)
            {
                result += entries.size();
            }
        }
        return result;
    }
private:
    void EnumHelper(Iter begin, Iter end)
    {
        if(begin == end)
            return;
        
        for(Iter next = begin;next != end; ++next)
        {
            substrings_.push_back(string(begin, next + 1));
        }
        EnumHelper(begin+1, end);
    }
    string input_;
    vector<string> substrings_;
};
int BinSq(int in)
{
    int left = 0, right = in +1;
    while(left != right)
    {
        auto mid = left + (right - left)/2;
        if (mid*mid > in)
        {
            right = mid;
        }
        else
        {
        left = mid +1;
        }
    }
    return left -1;
}
/*
 1,3,5,7
 2,4,6
 1,3,2,5,4,7,6
 */
int ZigZag(const vector<int>& input)
{
    std::vector<int> result;
    std::vector<int> zig(input.size(),1);
    std::vector<int> zag(input.size(),1);
    int max_result = 1;
    for(int i = 1; i != input.size();++i)
    {
        for(int j = 0; j!= i; ++j)
        {
            if(input[j] < input[i] && zag[j] +1 > zig[i])
            {
                zig[i] = zag[j] +1;
                if(zig[i] > max_result)
                    max_result = zig[i];
            }
            else if(input[j] > input[i] && zig[j] +1 > zag[i])
            {
                zag[i] = zig[j] +1;
                if(zag[i] > max_result)
                    max_result = zag[i];
            }
        }
    }
    
    return max_result;
}
int Fibonacci(int n)
{
    
    int fib_1 = 0, fib_2 = 1;
    int result = fib_1 + fib_2;
    for(int index = 2; index != n + 1;++index)
    {
        result = fib_1 + fib_2;
        fib_1 = fib_2;
        fib_2 = result;
    }
    return result;
}
struct TimeElement
   {
       __int64_t time_ = 0;
       int sensor_id_ = 0;
   };

class Chronicler
{
public:
       struct TimeElementComparator
    {
        bool operator()(const TimeElement& left, const TimeElement& right)
        {
            
            if(left.time_ > right.time_)
                return true;
            if(left.time_ == right.time_)
            {
                return left.sensor_id_ > right.sensor_id_;
            }
            return false;
        }
    };
    Chronicler(int sensors):cached_(sensors)
    {
        
    }
    void Normalize(const std::vector<TimeElement>& time_series)
    {
        for(const auto& time_element : time_series)
        {
            time_line_.push(time_element);
        }
        auto time_line = time_line_.top().time_;
        
        while(!time_line_.empty())
        {
            if(time_line == time_line_.top().time_)
            {
                cached_[time_line_.top().sensor_id_] = time_line_.top();
                time_line_.pop();
            }
            else
            {
                time_line = time_line_.top().time_;
                DrainQueue();
            }
                
        }
        
    }
    private:
    void DrainQueue()
    {
        
    }
    std::priority_queue<TimeElement,std::vector<TimeElement>,TimeElementComparator> time_line_{TimeElementComparator{}};
    std::vector<TimeElement> cached_;
};
class LongestPalindromic
{
public:
    LongestPalindromic(const string& input):input_(input)
    {
        for(auto first = input.cbegin();first != input.cend();++first)
        {
            auto last = input.cend();
            while(first != last )
            {
                auto palindrome_match = EndsWith(first, last, *first);
                
                if(palindrome_match != first)
                {
                    if(IsPalindrome(first+1, palindrome_match))
                    {
                        int thisLength = distance(first, palindrome_match+1);
                        result = std::max(thisLength,result);
                        if(first + thisLength == input_.cend())
                            return;
                        break;
                    }
                    
                }
                last = palindrome_match;
            }
            
        }
    }
    int GetResult()const
    {
        return result;
    }
private:
    bool IsPalindrome(string::const_iterator first, string::const_iterator last)
    {
        while(first != last && first != --last)
        {
            if(*first != *last)
                return false;
            ++first;
        }
        return true;
    }
    string::const_iterator EndsWith(string::const_iterator first, string::const_iterator last, char ch)
    {
        while(first != last && --last !=first)
        {
            if(*last == ch)
                return last;
        }
        return first;
    }
    string input_;
    int result = 1;
};
#include <stack>

class WaterMark
{
public:
    WaterMark(vector<int> bars):bars_(std::move(bars))
    {
      if(bars_.size()<3)
          return;
        cups_.push(0);
        for(auto next = 1; next != bars_.size();++next)
        {
            while(cups_.size() > 1 && bars_[next] > bars_[cups_.top()] )
            {
                cups_.pop();
            }
            if(cups_.size() == 1 &&  bars_[next] > bars_[cups_.top()] && cups_.top() + 1 == next)
            {
                cups_.pop();
            }
            cups_.push(next);
        }
    }
    int CalculateArea()
    {
        int result = 0;
        while(cups_.size()>1)
        {
            int rBound = cups_.top();
            cups_.pop();
            int lBound = cups_.top();
            int height = std::min(bars_[lBound],bars_[rBound]);
            ++lBound;
            for(int begin = lBound ; begin != rBound; ++begin)
            {
                if(bars_[begin] > height)
                {
                    height = bars_[begin];
                    continue;
                }
                result += (height - bars_[begin]);
            }
        }
        return result;
    }
private:
    stack<int> cups_;
    vector<int> bars_;
};
class HistoGram
{
public:
    void BuildStack()
    {
        if(bars_.size() ==0 )
            return;
        
        RectangleStart rect = {0,bars_[0]};
        rectangle_epochs_.push(rect);
        int max_area = 0;
        int next = 1;
       
        for(; next != bars_.size(); ++next)
        {
            
            while(rectangle_epochs_.size() && bars_[next] < rectangle_epochs_.top().height)
            {
                int width = next - rectangle_epochs_.top().position;
                int height = rectangle_epochs_.top().height;
                int area = width * height;
                if(area > max_area)
                {
                    max_area = area;
                    result.first = rectangle_epochs_.top().position;
                    result.second = next;
                }
                rectangle_epochs_.pop();
            }
            if(rectangle_epochs_.size() == 0)
            {
               rectangle_epochs_.push(RectangleStart{0,bars_[next]});
            }
            else if(bars_[next] > rectangle_epochs_.top().height)
            {
               rectangle_epochs_.push(RectangleStart{rectangle_epochs_.top().position + 1,bars_[next]});
            }
            
        }
        
        while(rectangle_epochs_.size())
        {
            int width = next - rectangle_epochs_.top().position;
            int height = rectangle_epochs_.top().height;
            int area = width * height;
            if(area > max_area)
            {
                max_area = area;
                result.first = rectangle_epochs_.top().position;
                result.second = next;
            }
            rectangle_epochs_.pop();
        }
    }
    HistoGram(vector<int> bars):bars_(std::move(bars))
    {
        if(bars_.size() ==0 )
            return;
        
        BuildStack();
    }
    std::pair<int, int> GetResult()
    {
        return result;
    }
private:
    struct RectangleStart
    {
        size_t position;
        int height;
    };
    stack<RectangleStart> rectangle_epochs_;
    vector<int> bars_;
   
    std::pair<int, int> result;
};
class ContiguousSum
{
public:
    ContiguousSum(vector<int> input):input_(move(input))
    {
        
    }
    bool CanSum(int target)
    {
        auto q_begin = input_.cbegin();
        int current_sum = 0;
        for(auto q_end = q_begin;q_end != input_.cend(); ++ q_end)
        {
            current_sum += *q_end;
            
            while (current_sum > target) {
                current_sum -= *q_begin;
                ++q_begin;
            }
            
            if(current_sum == target)
                return true;
        }
        
        return false;
    }
private:
    vector<int> input_;
};
int max_Steal(vector<int> numbers)
{
    if(numbers.size() ==0)
        return 0;
    if(numbers.size() ==1)
        return numbers[0];
    int prev_2 = numbers[0];
    int prev_1 = max(numbers[0],numbers[1]);
    
    for(int i = 2; i != numbers.size(); ++i)
    {
        if( prev_2 + numbers[i] > prev_1)
        {
            int temp = prev_1;
            prev_1 = prev_2 + numbers[i];
            prev_2 = temp;
        }
    }
        
    return max(prev_2,prev_1);
}

int pivotIndex(vector<int> numbers) {
    int left_sum = 0;
    int right_sum = 0;
    int sum = 0;
   for(int i = 0, j = numbers.size();i!=numbers.size();)
   {
       
       if(left_sum ==0 || left_sum < right_sum)
       {
           left_sum += numbers[i];
           ++i;
       }
       if(right_sum ==0 || right_sum < left_sum)
       {
           --j;
           right_sum += numbers[j];
       }
       if(left_sum == right_sum && i + 1 == j)
           return i ;
       
       if(i == numbers.size() || j == 0)
           break;
   }
   return -1;

    
}




size_t FindInflection(const vector<int>& input)
{
    size_t begin = 0;
    size_t end = input.size();
    while(begin != end)
    {
        auto mid = begin + (end - begin)/2;
        if(input[mid] > input[begin])
            begin = mid ;
        else
            end = mid;
    }
    return end;
    
}

class Event
{
public:
    enum class Type {Arrived, Departed};
public:
    Event(int year, Event::Type type):year_(year), type_(type)
    {
        
    }
    bool operator < (const Event& other)const
    {
        if(other.year_ == year_)
        {
            return type_ == Event::Type::Arrived && other.type_ == Event::Type::Departed  ;
        }
        return  year_ < other.year_;
    }
    bool IsArrival()const
    {
        return type_ == Event::Type::Arrived;
    }
private:
    int year_ = 0;
    Event::Type type_;
};
class NPR
{
public:
    NPR(const vector<Event>& registry):registry_(registry)
    {
        sort(registry_.begin(),registry_.end());
    }
    int MaxPopulation()
    {
        int result = 0;
        int current_count = 0;
        for(const auto& evt : registry_)
        {
            if(evt.IsArrival())
            {
                ++current_count;
                result = max(result,current_count);
            }
            else
            {
                --current_count;
            }
        }
        return result;
    }
private:
    vector<Event> registry_;
};
template <typename RandIter, typename BinaryPredicate>
RandIter Partition(RandIter first, RandIter last, BinaryPredicate pred)
{
    
    first = find_if_not(first, last, pred);
    if(first == last)
        return last;
    for(auto next = first+1; next!= last;++next)
    {
        if(pred(*next))
        {
            iter_swap(next, first);
            ++first;
        }
    }
    return first;
}
template <typename RandIter>
void QuickSort(RandIter first, RandIter last)
{
    if(first == last)
        return;
    auto middle = first + (last - first)/2;
    auto pivot = *middle;
    middle = Partition(first, last,[pivot](auto elem){return elem < pivot;} );
    auto middle2 = Partition(middle,last,[pivot](auto elem){return!( pivot < elem);});
    QuickSort(first, middle);
    QuickSort(middle2, last);
}
template <typename RandIter>
int QuickSelect(RandIter first, RandIter last, int rank)
{
    auto begin = first;
    while(first != last)
    {
        auto middle1 = first + (last - first)/2;
        auto pivot = *middle1;
        middle1 = Partition(first, last,[pivot](auto elem){return elem < pivot;} );
        auto middle2 = Partition(middle1,last,[pivot](auto elem){return!( pivot < elem);});
        auto rank_lower = middle1 - begin;
        if(rank < rank_lower)
        {
            last = middle1;
            continue;
        }
        auto rank_upper = middle2 - begin;
        if(rank < rank_upper)
            return *middle1;
        first = middle2;
    }
    throw invalid_argument{"Invalid input"};
}
class CloudWalk
{
public:
};
class TrafficVolume
{
public:
    TrafficVolume(int time, double volume):time_(time), volume_(volume)
    {
        
    }
    bool operator >= (const TrafficVolume& other)
    {
        if(*this == other)
            return true;
        
       if (volume_ == other.volume_)
           return time_ > other.time_;
        return volume_ > other.volume_;
    }
    bool operator == (const TrafficVolume& other)
    {
        return volume_ == other.volume_ && time_ == other.time_;
    }
    bool operator < (const TrafficVolume& other)
    {
        if (volume_ == other.volume_)
            return time_ < other.time_;
        
        return volume_ < other.volume_;
    }
private:
    int time_;
    double volume_;
    
};
/*
 1,2,3,4
 */
#include <list>
class SlidingWindowMax
{
public:
    SlidingWindowMax(vector<TrafficVolume> traffic ):traffic_(std::move(traffic))
    {
        
    }
    vector<TrafficVolume> GetMax(size_t window)
    {
        vector<TrafficVolume> result;
        for(auto first = traffic_.begin();first != traffic_.end();++first)
        {
            Push(first);
            while(distance(Max(), first) > window)
            {
                Pop();
            }
            result.push_back(*Max());
        }
        return result;
    }
private:
    using QIter = vector<TrafficVolume>::iterator;
    void Push(QIter iter)
    {
        while(!max_queue_.empty() && *max_queue_.back() < *iter)
        {
            max_queue_.pop_back();
        }
        max_queue_.push_back(iter);
    }
    void Pop()
    {
        max_queue_.pop_front();
    }
    QIter Max()
    {
        return *max_queue_.begin();
    }
    
    vector<TrafficVolume> traffic_;
    list<QIter> max_queue_;
};
class FindInWindow
{
public:
    FindInWindow(vector<string> text, vector<string> key_words):text_(move(text)), key_words_(move(key_words))
    {
        for(auto key_word : key_words_)
        {
            keyword_map_[key_word] = positions_list_.cend();
        }
    }
    pair<int,int> GetWindow()
    {
        pair<int,int> result;
        size_t minimum_window = std::numeric_limits<size_t>::max();
        for(auto first = text_.cbegin();first != text_.cend();++first)
        {
            if(keyword_map_.count(*first))
            {
               if(keyword_map_[*first] != positions_list_.cend())
               {
                   positions_list_.erase(keyword_map_[*first]);
               }
                positions_list_.push_back(first);
                keyword_map_[*first] = --positions_list_.cend();
                if(positions_list_.size() == key_words_.size())
                {
                    size_t this_distance = *(--positions_list_.cend()) - *positions_list_.cbegin();
                    if(this_distance < minimum_window)
                    {
                        minimum_window = this_distance;
                        result = {distance(text_.cbegin(),*positions_list_.cbegin()),
                            distance(text_.cbegin(),*(--positions_list_.cend()))
                        };
                    }
                }
            }
        }
        return result;
    }
private:
    vector<string> text_;
    vector<string> key_words_;
    list<vector<string>::const_iterator> positions_list_;
    map<string,list<vector<string>::const_iterator>::const_iterator> keyword_map_;
    
};
bool IsVowel(char ch)
{
    return ch == 'a' || ch == 'e' || ch == 'i' || ch == 'o' || ch == 'u';
}
void ReverseVowels(string& input)
{
    auto begin = input.begin();
    auto end = input.end();
    while(begin != end)
    {
        if(IsVowel(*begin) && IsVowel(*(end -1)))
        {
            //swap(input[begin],input[end-1]);
            iter_swap(begin++,--end);
        }
        if(!IsVowel(*begin))
        {
            ++begin;
        }
        if(!IsVowel(*(end -1)))
        {
            --end;
        }
    }
}
bool CanPalindrome(const string& input)
{
    
    bool result = false;
    map<char,int> freq;
    for(const auto& ch : input)
    {
        ++freq[ch];
    }
    int count_of_1 = 0;
    int count_of_2 = 0;
    int count_of_3 = 0;
    
    for (const auto& elem:freq )
    {
        if(elem.second == 1)
        {
            ++count_of_1;
        }
        else if(elem.second == 2)
        {
            ++count_of_2;
        }
        else if(elem.second == 3)
        {
            ++count_of_3;
        }
    }

    if(input.size() % 2)
    {
        if(count_of_2*2 + count_of_1 == input.size())
            return true;
        
        if(count_of_2*2 + 3*count_of_3 == input.size())
            return true;
        
        return false;
        //odd so can be made even length pali
        // abba, and one char can have 3, or one, size -1 chars to have 2 counts
        
    }
    else
    {
        //even so can be made odd length
        //liril -> 2 chars with one count, 1 char with two counts or one char with 3 counts and rest with 2 counts
        
        // all with 2 counts
        // 2 chars have 1 count, rest all 2 counts
        //1 char 3, rest all 2
        //2 chars with one count && rest all with 2 counts
        //1 char with 3 counts, 1 char 1 count, rest with 2 counts
        //1 char with 2 counts || 1 char with 3 counts, && size -1 chars with 2 counts
        if(count_of_2*2 == input.size())
            return true;
        
        if(count_of_2*2 + count_of_1 + 3*count_of_3 == input.size())
            return true;
        if(count_of_2*2 + 2*count_of_1  == input.size())
            return (count_of_1);
    }
           
    return result;
}
#include <array>
class GridLocked
{
public:
    GridLocked(vector<vector<int>> grid):grid_(move(grid))
    {
        i_max_ = grid_.size();
        j_max_ = grid_[0].size();
    }
    void FillGrid()
    {
        for(int i = 0; i != i_max_; ++i)
            for(int j =0; j!= j_max_; ++j)
                if(IsRoom(i, j))
                    DoBFS(i, j);
    }
private:
    void DoBFS(int i, int j)
    {
        vector<int> distance(i_max_*j_max_,-1);
        queue<pair<int, int>> bfs_queue;
        bfs_queue.push({i,j});
        distance[i*j_max_ + j] = 0;
        while(!bfs_queue.empty())
        {
        
            auto src_elem = bfs_queue.front();
            bfs_queue.pop();
            auto neighbors = GetNeighbors(src_elem.first, src_elem.second);
            for(const auto& neighbor: neighbors)
            {
                if(IsGate(neighbor.first, neighbor.second))
                {
                    
                    grid_[i][j] = distance[src_elem.first*j_max_ + src_elem.second] + 1;
                    return;
                }
               if((distance[neighbor.first*j_max_ + neighbor.second] == -1) && !IsWall(neighbor.first, neighbor.second))
                {
                    distance[neighbor.first*j_max_ + neighbor.second] = distance[src_elem.first*j_max_ + src_elem.second] + 1;
                    bfs_queue.push({neighbor.first,neighbor.second});
                }
            }
        }
        
    }
    vector<vector<int>> grid_;
    size_t i_max_ = 0;
    size_t j_max_ = 0;
    bool IsWall(int i, int j)const
    {
        return grid_[i][j] == CellType::Wall;
    }
    bool IsRoom(int i, int j)const
    {
        return grid_[i][j] == CellType::Room;
    }
    bool IsGate(int i, int j)const
    {
        return grid_[i][j] == CellType::Gate;
    }
    vector<pair<int,int>> GetNeighbors(int i, int j)
    {
        vector<pair<int,int>> result;
        for(auto offset : offsets_)
        {
            if(IsValidIrange(offset.first + i) && IsValidJrange(offset.second +j))
            {
                result.push_back({offset.first + i,offset.second +j});
            }
        }
        return result;
    }
    array<pair<int, int>,4> offsets_ = {
        {
        {-1,0},
        {1,0},
        {0,-1},
        {0,1}
        }
    };
    bool IsValidIrange(int i) const
    {
        if(i < 0)
            return  false;
        if(! (i < i_max_  ))
            return false;
        return true;
    }
    bool IsValidJrange(int j) const
    {
        if(j < 0)
            return  false;
        if(! (j < j_max_  ))
            return false;
        return true;
    }
    enum CellType {Wall = -1, Room = std::numeric_limits<int>::max(), Gate = 0};
    
};

int MaxBoats( vector<int>& riders, int limit)
{
    int result = 0;
    sort(riders.begin(),riders.end());
    auto being = riders.begin();
    auto end = riders.end();
    auto count = riders.size();
    while (being != end)
    {
        if(limit < *being + *(end -1 ))
        {
            ++result;
            --end;
            --count;
        }
        else
        {
            ++result;
            ++being;
            --end;
            count -=2;
        }
    }
    if(count)
        ++result;
    return result;
}
class MaxPartition
{
public:
    MaxPartition(string input):input_(move(input))
    {
        
    }
    vector<size_t> PartitionUnique()
    {
        vector<size_t> result;
        map<char, pair<size_t,size_t>> positions_map;
        size_t pos = 0;
        size_t end = input_.size();
        for(auto ch : input_)
        {
            if(positions_map.count(ch) == 0)
            {
                positions_map[ch] = {pos,pos};
            }
            else
            {
                positions_map[ch].second = pos;
            }
            ++pos;
        }
        size_t cursor = 0;
        
        while(cursor != end)
        {
            //ababcbacadefegdehijhklij
            {
                auto partition_begin = cursor; //positions_map[input_[cursor]].first;
                auto partition_end = positions_map[input_[cursor]].second;
                for(auto next = partition_begin + 1;next != partition_end; ++next)
                {
                    if(positions_map[input_[next]].second > partition_end)
                        partition_end = positions_map[input_[next]].second;
                }
                result.push_back(partition_end - partition_begin + 1);
                cursor = partition_end + 1;
            }
        }
        return  result;
    }
private:
    string input_;
};

/*
 100
 
 */


string HextToString(int input)
{
    string result;
    while(input)
    {
        auto ch = input % 16;
        if(ch > 9)
        {
            ch = ch - 10  + 'A';
        }
        else
        {
            ch = ch + '0';
        }
        result.push_back(ch);
        input /= 16;
    }
    reverse(result.begin(), result.end());
    return result;
}
string NextTime(string input)
{
    string result;
    if(input[3] == '5' && input[4] == '9')
    {
        auto ch = min(input[0],input[1]);
        result += ch;
        result +=ch;
        result += ':';
        result += ch;
        result +=ch;
        
    }
    else
    {
        //if(
        auto ch = min(input[0],input[1]);
        if(ch > input[4])
        {
            result = input;
            result[4] = ch;
            return result;
        }
        auto ch_max = max(input[0],input[1]);
        if(ch_max > input[4])
        {
            result = input;
            result[4] = ch_max;
            return result;
        }
        
    }
    return  result;
}
bool AnagramPalindrome(string word) {
  
    //lirill
  map<char, int> word_freq;
  for(auto ch : word)
  {
      ++word_freq[ch];
  }
  int count_of_odds = 0;
  int count_of_evens = 0;
  for(auto elem : word_freq)
  {
     if(elem.second % 2)
     {
         ++count_of_odds;
     }
     else if(elem.second == 2)
     {
         ++count_of_evens;
     }
  }
    if((count_of_odds == 2) || (count_of_odds ==1))
        return true;
    return false;
}

int findFulcrum(vector<int> numbers) {
    
    vector<int> partial_sum;
    //1,2,3
    //0,1,3,6
    int sum_local = 0;
    partial_sum.push_back(sum_local);
    for(auto elem: numbers)
    {
        sum_local += elem;
        partial_sum.push_back(sum_local);
    }
    for(int index = 1; index != partial_sum.size(); ++index)
    {
        int left_sum = partial_sum[index -1 ];
        int right_sum = sum_local - left_sum - numbers[index-1];
        if(left_sum == right_sum)
            return index -1;
    }
    return -1;
}
class PowerSet
{
public:
    PowerSet(string input):input_(input)
    {
        
    }
    vector<string> ComputePowerSet()
    {
        string temp_result;
        PowerSetHelper(input_.begin(),input_.end(),temp_result);
        return result;
    }
private:
    void PowerSetHelper(string::iterator cursor, string::iterator end, string& this_set)
    {
      if(cursor == end)
      {
          result.push_back(this_set);
          return;
      }
      
      this_set.push_back(*cursor);
      PowerSetHelper(cursor + 1, end, this_set);
      this_set.pop_back();
        
        
    }
    string input_;
    vector<string> result;
};
class LIS
{
public:
    LIS(vector<int> input):input_(move(input))
    {
        input_length_ = input_.size();
    }
    size_t GetResult()
    {
        if(input_length_ <2)
        {
            return input_length_;
        }
        vector<size_t> dp(input_length_,1);
        size_t max_so_far = std::numeric_limits<size_t>::min();
        for(size_t index = 1; index != input_length_; ++index)
        {
            for(size_t pred = 0; pred != index; ++pred)
            {
                if(input_[pred] < input_[index])
                {
                    if(dp[index] < dp[pred] + 1)
                    {
                        dp[index] = dp[pred] + 1;
                        if(max_so_far < dp[index])
                        {
                             max_so_far = dp[index];
                        }
                    }
                }
            }
        }
        return max_so_far;
        
    }
private:
    vector<int> input_;
    size_t input_length_;
};
// 3[ab2[a]]abcd
class LCS
{
public:
    LCS(string left, string right):left_(left), right_(right)
    {
        
    }
    size_t Compute()
    {
        vector<size_t> dp(right_.length() + 1);
        for(size_t rhs = 1; rhs != right_.length() + 1;++rhs)
        {
            size_t left_top = 0;
            for(size_t lhs =0; lhs != left_.length(); ++lhs)
            {
             size_t next_left_top= next_left_top = dp[rhs];
             if(left_[lhs] != right_[rhs-1])
             {
                 auto max_new = max(left_top,dp[rhs-1]);
                 dp[rhs] = max(max_new,dp[rhs]);
             }
            else
            {
                dp[rhs] = left_top + 1;
            }
                left_top = next_left_top;
            }
        }
        return dp[right_.length()];
    }
private:
    string left_;
    string right_;
    
};
class StockMaster
{
public:
    StockMaster(vector<size_t> stocks,int k ):stocks_(stocks),k_(k)
    {
        
    }
    //5,11,3,50,60,90
    vector<size_t> GetTransactions()
    {
        vector<size_t> result;
        size_t min_cursor = stocks_[0];
        
        size_t max_profit = 0;
        for(auto stock: stocks_)
        {
            if(stock < min_cursor )
            {
                if(max_profit)
                {
                    result.push_back(max_profit);
                    max_profit = 0;
                }
                min_cursor = stock;
            }
            else if(stock - min_cursor > max_profit)
            {
                max_profit = stock - min_cursor;
            }
                
        }
        if(max_profit)
            result.push_back(max_profit);
        return result;
    }
    size_t MaxProfitWithCooldown()
    {
        
        vector<size_t> dp(stocks_.size(),0);
        for(int i = 1; i != stocks_.size(); ++i)
        {
            dp[i] = dp[i-1];
            for(int j = 0; j !=i; ++j)
            {
            //1. can sell ?
                if(stocks_[j]< stocks_[i])
                {
                    if(j < 2 == false)
                    {
                    
                        dp[i] = max(dp[i],dp[j-2] + stocks_[i] - stocks_[j]);
                    }
                    else
                    {
                        dp[i] = max(dp[i],stocks_[i] - stocks_[j]);
                    }
                }
            }
        }
        return dp.back();
    }
    size_t MaxProfit()
    {
        vector<size_t> dp(stocks_.size());
        for(int i = 1; i != stocks_.size(); ++i)
        {
            dp[i] = dp[i-1];
            for(int j =0; j!= i ;++j)
            {
                if(stocks_[j] < stocks_[i])
                {
                    dp[i] = max(dp[i],(j ? dp[j-1] :0)+ stocks_[i]- stocks_[j]);
                }
            }
        }
        return dp.back();
    }
private:
    vector<size_t> stocks_;
    int k_;
};

class TreeNode
{
public:
    TreeNode(int data,shared_ptr<TreeNode> left, shared_ptr<TreeNode>right):
    left_{left},right_{right}
    {
        
    }
    void AddLeft(shared_ptr<TreeNode> left)
    {
        left_ = left;
    }
    void AddRight(shared_ptr<TreeNode> right)
    {
        right_ = right;
    }
    shared_ptr<TreeNode>& Left()
    {
        return left_;
    }
    shared_ptr<TreeNode>& Right()
    {
        return right_;
    }
    int Data()
    {
        return  data_;
    }
private:
    int data_;
    shared_ptr<TreeNode> left_,right_;
    
};
class LNode
{
    
};
vector<shared_ptr<TreeNode>> EnumTrees(int count)
{
    
    if(!count)
    {
        vector<shared_ptr<TreeNode>> result{nullptr};
        return result;
    }
    vector<shared_ptr<TreeNode>> result;
    for(auto begin = 0; begin != count; ++begin)
    {
        auto left_count = begin - 0;
        auto right_count = count - (begin + 1);
        auto left_trees = EnumTrees(left_count);
        auto right_trees = EnumTrees(right_count);
        for(auto left_tree :left_trees)
            for(auto right_tree : right_trees)
            {
                result.push_back(make_shared<TreeNode>(0,move(left_tree),move(right_tree)));
            }
    }
    return  result;
    
}
shared_ptr<TreeNode> MakeTree(shared_ptr<TreeNode>& head,int begin, int end)
{
    if(begin == end)
        return nullptr;
    auto mid = begin + (end - begin)/2;
    auto left_subtree = MakeTree(head, begin, mid);
    auto root = head;
    root->Left() = left_subtree;
    head = head->Right();
    root->Right() = MakeTree(head, mid+1, end);
    return root;
}
shared_ptr<TreeNode> MakeTree(shared_ptr<TreeNode> head,int count)
{
    return MakeTree(head,0, count);
}
//ADOBECODEBANC, ABC
struct Matched
    {
     int expected = 0;
     int found = 0;
    };

class MinWindow
{
public:
MinWindow(string input, string pattern):input_(input),pattern_(pattern)
{
    for(auto ch: pattern_)
    {
        ++matches[ch].expected;
    }
    target_count = pattern_.length();
}
/*
  A A B C D D
A 1 1 1 1 1 1
C 1 1 1 2 2 2
D 1 1 1 2 3 3
 */
string GetMin()
{
    string result;
    
    
    int left = 0;
    
    int result_left = 0;
    int result_right = 0;
    int min_so_far = numeric_limits<int>::max();
    left = Begin();
    if(left == input_.length())
        return "";

    for(auto index = left; index != input_.length();++index)
    {
        if(matches.count(input_[index]))
        {
            ++matches[input_[index]].found;
            if(matches[input_[index]].found == matches[input_[index]].expected)
            {
                  --target_count;
                  if(target_count ==0)
                  {
                      if(index + 1 - left < min_so_far )
                      {
                          result_left = left;
                          result_right = index + 1;
                          min_so_far = result_right -  result_left;
                      }
                      left = Drain(left);
                  }
            }
        }
        
    }
    //while(index - )
    return string(input_.begin() + result_left,input_.begin() + result_right);
}
private:
int Drain(int left)
{
  --matches[input_[left]].found;
    ++target_count;
  for(++left;left != input_.length();++left)
  {
      if(matches.count(input_[left]))
         return left;
  }
  return left;
}
int Begin()
{
    for(auto index = 0; index != input_.length();++index)
    {
        if(matches.count(input_[index]))
         return index;
    }
    return input_.length();
}
string input_;
string pattern_;
map<char,Matched> matches;
    int target_count;
};
string minimumWindowSubstring(string fullString, string chars) {
MinWindow mw{fullString,chars};
return mw.GetMin();
}
/*
    _ K D B C C
 _  1 1 1 1 1 1
 A  0 1 1 1 1 1
 B  0 0 0 1 1 1
 C  0 0 0 0 1 2
 */
int CountDistinct(string input, string text)
{
    /*
     if(src[i] == input[j]
     {
       dp[i] = dp[i-1][j-1] + dp[i][j-1]
     }
     else
     {
      dp[i] = dp[i][j-1]
     }
     */
    vector<vector<int>> dp(input.length() + 1,vector<int>(text.length() + 1));
    dp[0] = vector<int>(vector<int>(text.length() + 1,1));
    
    for(int i = 1; i != input.length() + 1; ++i)
    {
        for(int j = 1; j != text.length() + 1; ++j)
        {
            if(input[i-1] == text[j-1])
            {
                dp[i][j] = dp[i-1][j-1] + dp[i][j-1];
            }
            else
            {
                dp[i][j] = dp[i][j-1];
            }
        }
    }
    
    return dp.back().back();
    
}

//
class Paranthesize
{
public:
    Paranthesize(int count)
    {
        GenerateAll(count, count);
    }
    
private:
    void GenerateAll(int left_bal , int right_bal)
    {
        if(right_bal == 0 && left_bal == 0)
        {
            results_.push_back(result_);
            return;
        }
        if(left_bal)
        {
            result_.push_back('(');
            GenerateAll(left_bal -1, right_bal);
            result_.pop_back();
        }
        if(right_bal > left_bal)
        {
            result_.push_back(')');
            GenerateAll(left_bal, right_bal - 1 );
            result_.pop_back();
        }
    }
    string result_;
    vector<string> results_;
};
class PalindromicDecompositon
{
public:
    PalindromicDecompositon(string input):input_(input)
    {
        last_ = input_.cend();
        GenerateRecursively(input_.cbegin());
    }
private:
    bool  GenerateRecursively(string::const_iterator cursor )
    {
        if(cursor == last_)
            return true;
        
        for(auto next = cursor;next != last_; ++next)
        {
            if(IsPalindrome(cursor, next + 1))
            {
                if(GenerateRecursively(next + 1))
                {
                    result_.push_back(string{cursor,next +1});
                }
            }
            else if(next + 1 == last_)
            {
                return false;
            }
        }
        return true;
    }
    
    bool IsPalindrome(string::const_iterator first, string::const_iterator last)
    {
        while(first != last && first != --last)
        {
            if(*first != *last)
                return false;
            ++first;
        }
        return true;
    }
    vector<string> result_;
    string input_;
    string::const_iterator last_;
};
class Solution {
public:
    size_t uniquePathsWithObstacles(vector<vector<int>>& obstacleGrid) {
        
        
        size_t y = obstacleGrid.size();
        size_t x = obstacleGrid[0].size();
        vector<vector<size_t>> dp(y,vector<size_t>(x));
        for(int i = 0; i != y; ++i)
        {
            for(int j = 0; j != x; ++j)
            {
                dp[i][j] = !obstacleGrid[i][j];
            }
        }
        for(int j = 1; j != x; ++j)
        {
            if(dp[0][j])
            {
                dp[0][j] = dp[0][j -1];
            }
        }
        for(int i = 1; i != y; ++i)
        {
            if(dp[i][0])
            {
                dp[i][0] = dp[i-1][0];
            }
        }
        for(int i = 1; i != y; ++i)
        {
            for(int j = 1;  j != x; ++j)
            {
                if(dp[i][j])
                {
                    dp[i][j] = dp[i-1][j] + dp[i][j-1];
                }
            }
        }
        
        return dp.back().back();
    }
};
/*
 
  0 1 2 3 4  5
1 0 0 1 2 3  4
2 0 0 0 0 1  2
3            0
4            0
 */
#include <numeric>
class SplitFairly
{
public:
    SplitFairly(vector<int> input):input_(input)
    {
        target_ = accumulate(input_.cbegin(),input_.cend(),0);
    }

    vector<int> Compute()
    {
        vector<vector<int>> dp(input_.size(),vector<int>(target_/2 + 1) );
        vector<int> result;
        for(int capacity = target_/2;  capacity !=0 ;--capacity)
        {
            if(capacity < input_[0] == false)
                dp[0][capacity] = capacity - input_[0];
            else
            {
                dp[0][capacity] = capacity;
            }
        }
        for(int item = 1; item != input_.size(); ++ item)
        {
            for(int capacity = target_/2; capacity != 0  ;--capacity)
            {
                int wo_item = dp[item-1][capacity];
                if(capacity < input_[item] == false)
                {
                    int with_item = dp[item-1][capacity - input_[item]];
                    dp[item][capacity] = min(with_item,wo_item);
                }
                else
                {
                    dp[item][capacity] = wo_item;
                }
            }
        }
        vector<int> result1;
        {
            
            int j = target_/2;
            int i = input_.size()-1;
            
            while(i < 0 == false)
            {
                if(i ==0 || dp[i][j] != dp[i-1][j])
                {
                    if(j < input_[i] == false)
                    {
                        result.push_back(input_[i]);
                        j -= input_[i];
                    }
                    
                }
                else
                {
                    result1.push_back(input_[i]);
                }
                --i;
            }
         
        }
        return result;
    }
private:
    int target_;
    vector<int> input_;
};
class RegExp
{
public:
    RegExp(string pattern):pattern_(pattern)
    {
        
    }
    bool Match(string text)
    {
        vector<vector<int>> dp(text.size() + 1,vector<int>(pattern_.size() + 1));
        dp[0][0] = 1;
        //*
        for(int j = 1; j !=pattern_.size() + 1; ++j )
        {
            if(pattern_[j-1] == '*')
            {
                dp[0][j] = dp[0][j-1] || dp[0][j-2];
            }
          //dp[0][]
        }
        for(int i = 1; i != text.size() + 1; ++i )
        {
            dp[i][0] = 1;
            for(int j = 1; j != pattern_.size() + 1; ++j)
            {
                if(text[i-1] == pattern_[j-1])
                {
                    dp[i][j] = dp[i-1][j-1];
                }
                else
                {
                    //ab*
                    //.*
                    if(pattern_[j-1] == '*')
                    {
                        dp[i][j] = (j -1 ==0)||(j -2 == 0) || dp[i-1][j] ||
                        dp[i][j-2] || pattern_[j-2] == '.' || pattern_[j-2] == text[i-1];
                    }
                }
            }
        }
        return dp.back().back() == 1;
    }
private:
    string pattern_;
};
long repeatedString(string s, long n) {

    long result = 0;
    long fullStrings = n/s.size();
    long count = 0;
    for(int i =0; i!= s.size(); ++i)
    {
        if(s[i] == 'a')
            ++count;
    }
    long remainder = n % s.size();
    for(size_t i = 0; i != remainder; ++i)
    {
        if(s[i] == 'a')
           ++result;
    }
    
    return  fullStrings*count +  result;
}
int numberPairs(int n, vector<int> ar) {

    map<int,int> histogram;
    for(auto elem : ar)
    {
        ++histogram[elem];
    }
    int result = 0;
    for( auto& elem : histogram)
    {
        result += elem.second/2;
    }
    return result;
}
template<typename T>
class ConcurrentVector
{
public:
    ConcurrentVector()
    {
        
    }
    T& operator[](size_t index)
    {
        TRef ref{this,index};
        return move(ref);
    }
    void push_back(const T& val)
    {
        vector_.push_back(val);
    }
    class TRef
    {
    public:
        TRef(ConcurrentVector<T>* owner, size_t index):owner_(owner),sc_lk_(owner->mtx_),
        index_(index)
        {
            
        }
        operator T&()
        {
            return owner_->vector_[index_];
        }
        ~TRef()
        {
            owner_ = nullptr;
        }
        TRef(const TRef&) = delete;
        TRef( TRef&&) = default;
    private:
        ConcurrentVector<T> *owner_;
        scoped_lock<mutex> sc_lk_;
        size_t index_;
    };
private:
    vector<int> vector_;
    mutex mtx_;
    
};
void AFunc(int& input)
{
    ++input;
}
/*
  _ b b
_ 1 0 0
b 1 1 0
b 1 2 1
b 1 3 3
 */
class CountSubsequences
{
public:
    CountSubsequences(string src, string dest):src_(src), dest_(dest)
    {
        
    }
    size_t Count()
    {
        vector<size_t> dp(dest_.size());
        //dp[0] = 1;
        size_t diagonal = 1;
        for(size_t i = 0; i != src_.size() ; ++ i)
        {
            for(size_t j = 0; j != dest_.size() ; ++j)
            {
                size_t next = dp[j];
                if(src_[i] == dest_[j])
                {
                    dp[j] += diagonal;
                }
                diagonal = next;
            }
            diagonal = 1;
            
        }
        return dp.back();
    }
private:
    string src_;
    string dest_;
};
int countingValleys(int n, string s)
{
    int step_sum = 0;
    int result = 0;
    for(auto ch : s)
    {
        if(ch == 'U')
        {
            ++step_sum;
            if(step_sum == 0)
            {
                ++result;
            }
        }
        else
        {
            --step_sum;
        }
    }
    
    return result;
}
//0,1,0,0,0,1,0
int jumpingOnClouds(vector<int> c) {

    int steps = 0;
    for(size_t i = 0; i != c.size(); )
    {
        if(c[i] ==0)
        {
            if(i + 1 == c.size() - 1 || i + 2 == c.size()-1)
            {
                return ++steps;
            }
            i += c[i + 2] == 0 ? 2 : 1;
            ++steps;
        }
    }
    return 0;
}
void checkMagazine(vector<string> magazine, vector<string> note)
{
    map<string,int> note_histogram;
    for(auto note_word : note)
    {
        ++note_histogram[note_word];
    }
    for(auto mag_word : magazine)
    {
        if(note_histogram.count(mag_word))
        {
            if(--note_histogram[mag_word] ==0)
            {
                note_histogram.erase(mag_word);
                if(note_histogram.size() ==0)
                {
                    cout << "Yes";
                    return;
                }
            }
        }
    }
    cout << "No";
}
//cde, abc -> cde and abc

size_t makeAnagramEx(string a, string b) {
     
    sort(a.begin(),a.end());
    sort(b.begin(),b.end());
    auto first1 = a.begin();
    auto last1 = a.end();
    
    auto first2 = b.begin();
    auto last2 = b.end();
    size_t count = 0;
    
    while(first1 != last1)
    {
        if(first2 == last2)
        {
            break;
        }
        if(*first2 < *first1)
        {
            ++count;
            ++first2;
        }
        else if(*first2 == *first1)
        {
            ++first1;
            ++first2;
        }
        else
        {
            ++first1;
            ++count;
        }
    }
    return count + distance(first1, last1) + distance(first2, last2);
}
int hourglassSum(vector<vector<int>> arr)
{
    if(arr.size() < 3 || arr[0].size() < 3)
        return 0;
    
    int result = 0;
    for(int i = 1; i != arr.size() -1 ; ++i)
    {
        for(int j = 1; j != arr[0].size() -1; ++j)
        {
            auto this_hour_glass_sum = arr[i][j];
            this_hour_glass_sum += arr[i-1][j-1];
            this_hour_glass_sum += arr[i-1][j];
            this_hour_glass_sum += arr[i-1][j+1];
            
            this_hour_glass_sum += arr[i+1][j-1];
            this_hour_glass_sum += arr[i+1][j];
            this_hour_glass_sum += arr[i+1][j+1];
            result = max(result,this_hour_glass_sum);
        }
        //for(in)
    }
    return result;
}
long minTime(vector<long> machines, long goal)
{
    long days = 1;
    long produced = 0;
    while(produced < goal)
    {
        produced = 0;
        for(auto days_per_item : machines)
        {
            produced += days / days_per_item;
        }
    }
    return days;
}
long largestRectangle(vector<int> h) {
    
    long max_area = 0;
    stack<int> dp;
    for(int i = 0; i != h.size(); ++i)
    {
        while(!dp.empty() &&  h[i] > h[dp.top()] == false)
        {
            long local_area = h[dp.top()] * (i - dp.top());
            max_area = max(local_area,max_area);
            dp.pop();
        }
        dp.push(i);
    }
    //
    
    auto right = h.size();
    while(!dp.empty())
    {
        auto ht = h[dp.top()];
        dp.pop();
        auto left = dp.empty() ? 0 : dp.top() + 1;
        long local_area = ht *(right - left);
        max_area = max(local_area,max_area);
    }
    return max_area;
}
string isBalanced(string s) {
    
    stack<char> brace;

    for(auto ch : s)
    {
        switch (ch) {
            case '{':
            case '[':
            case '(':
                brace.push(ch);
                break;
                
            case '}':
                if(brace.empty() || brace.top() != '{')
                    return "NO";
                brace.pop();
                
                break;
                
            case ']':
                if(brace.empty() || brace.top() != '[')
                    return "NO";
                brace.pop();
                
                break;
                
            case ')':
                if(brace.empty() || brace.top() != '(')
                    return "NO";
                brace.pop();
                break;
            default:
                break;
        }
    }
    return brace.empty() ? "YES" : "NO";
}
int DoBFS(vector<vector<int>>& grid, int i, int j)
{
    int result = 1;
    queue<pair<int,int>> bfs;
    bfs.push({i,j});
    while(!bfs.empty())
    {
        int left = j > 0 ? grid[i][j-1]:0;
        //int top =
    }
    return 0;
}
int maxRegion(vector<vector<int>> grid) {
//queue<pa
  
for(int i = 0; i != grid.size(); ++ i)
    for(int j = 0; j != grid[0].size(); ++j)
    {
       if(grid[i][j] == 1)
       {
           grid[i][j] = 0;
           int count = DoBFS(grid,i,j);
       }
    }
    return 0;
}
class AnaGramm
{
public:
    AnaGramm(string text, string pattern):text_(move(text)), pattern_(pattern),pattern_size_(pattern.size())
    {
        for(auto p: pattern)
        {
            ++histograph_[p];
        }
        
    }
    size_t Find()
    {
        
        for(size_t i = pattern_size_; i < text_.size() + 1;)
        {
            --i;
            if(histograph_.count(text_[i]))
            {
                if(FindMatch(i + 1))
                return i - pattern_size_;
                i+=2;
            }
            else
            {
                ++i;
                i +=pattern_size_;
            }
        }
        return text_.size();
    }
    bool FindMatch(size_t i)
    {
        auto begin = i - pattern_size_;
        while( i != begin)
        {
            --i;
            if(histograph_.count(text_[i]) &&
               ( text_histograph_[text_[i]]  < histograph_[text_[i]] ))
            {
                ++text_histograph_[text_[i]];
            }
            else
            {
                for(i = i + 1;i != begin + pattern_size_; ++i )
                {
                    --text_histograph_[text_[i]];
                }
                return false;
            }
        }
        return true;
    }
private:
    string text_, pattern_;
    const size_t pattern_size_;
    map<char,int> histograph_;
    map<char,int> text_histograph_;
};
class Grid
{
    
};
class OptimumTravel
{
public:
    enum Ordering {SMALLER,EQUAL,LARGER};
    Ordering Compare(double x, double y)
    {
        auto diff = (x - y)/y;
        return diff < - numeric_limits<double>::epsilon() ? SMALLER
        : diff < numeric_limits<double>::epsilon() ? LARGER : EQUAL;
    }
    
    OptimumTravel(vector<int> width, vector<int> speed,int offset )
    :width_{width},speed_{speed},offset_{offset}
    {
        
    }
    
    double GetTime(double x, double y, int speed)
    {
        auto dist = std::sqrt(x*x + y*y);
        return dist/speed;
    }
    
    double GetSlope(double left, double right, size_t index)
    {
        auto dx = right - left;
        auto dy = 0.0;
        for(auto i = index; i != width_.size(); ++i)
        {
            dy += width_[i];
        }
        auto slope = dy/dx;
        return slope;
    }
    double GetDiagonalTime(double left, double right, size_t index)
    {
        auto slope = GetSlope(left,right,index);
        double time = 0.0;
        for(auto i = index; i != width_.size(); ++i)
        {
            auto x = width_[i]/slope;
            auto y = width_[i];
            auto dist = std::sqrt(x*x + y*y);
            time += dist/speed_[i];
        }
        return time;
    }
    double ComputeTime(size_t begin, size_t end,int left, int right)
    {
        if(begin + 1 == end)
        {
            double y = width_[begin];
            double x = right - left;
            return GetTime(x,y,speed_[begin]);
        }
        auto ref_time = GetDiagonalTime(left,right,begin);
        auto range_begin = left;
        auto range_end = right;
         
        while(range_begin < range_end)
        {
            auto mid = range_begin + (range_end - range_begin)/2;
            auto my_time = GetTime(mid - left, width_[begin], speed_[begin]);
            double top_strips_time = ComputeTime(begin + 1, end, mid, right);
            auto total_time = my_time + top_strips_time;
            
            if(Compare(total_time, ref_time) == Ordering::LARGER)
            {
                range_end = mid;
            }
            else if(Compare(total_time, ref_time) == Ordering::SMALLER)
            {
                ref_time = total_time;
                range_begin = mid;
            }
            else
                return ref_time;
        }
        return ref_time;
    }
    double ComputeTime()
    {
        auto result = ComputeTime(0,width_.size(),0,offset_);
        return result;
    }
    
private:
    vector<int> width_;
    vector<int> speed_;
    int offset_;
    vector<double> delta_x_;
};
class MaxPath
{
public:
    MaxPath(shared_ptr<TreeNode> root):root_(root)
    {
        
    }
int ComputePath(shared_ptr<TreeNode> root)
{
    if(root == nullptr)
    {
        return 0;
    }
    auto left = ComputePath(root->Left());
    auto right = ComputePath(root->Right());
    auto sum_single_leg = max(max(left,right) + root->Data(),root->Data());
    max_global = max(max_global,max(sum_single_leg,left + right + root->Data()));
    return sum_single_leg;
}
private:
    int max_global = 0;
    
    shared_ptr<TreeNode> root_;
};
/*
 #include <iostream>
 using namespace std;

 // hello

 struct Node
 {
     
     int data ;
     Node* next_ptr_ = nullptr;
 };
 
 //pre condition - list_head is non null ,
 bool RemoveNode(Node** list_head_ptr,Node* node_to_remove)
 {
   assert(list_head);
   assert(node_to_remove);
   assert(list_head_ptr);
   Node* list_head = *list_head_ptr;
   Node* prev = nullptr;
   while(list_head && list_head != node_to_remove)
   {
     prev = list_head;
     list_head = list_head->next_ptr;
   }
   if(list_head == nullptr)
   {
     //failed to find the node_to_remove;
     //assert(list_head);
     return false;
   }
   if(prev)
   {
     prev->next =   prev->next->next;
   }
   else if(list_head->next)
   {
     
     *list_head_ptr =  list_head->next;
   }
   else
   {
     *list_head_ptr = nullptr;
   }
   return true;
 }

 int main() {
  Node * list = MakeSequentialList(5); // 1->2->3->4->5->NULL
     
     Node * oneNode = FindNode(list, 1);
    
     RemoveNode(list, oneNode);
     
     
     /// list: 2->3->4->5->NULL
     
   return 0;
 }

 */
class WorkBalancer
{
public:
    WorkBalancer(vector<int> work,int workers):
    work_(work),
    workers_(workers)
    {
        
    }
    int MinimizeLoad()
    {
        auto min_load = *std::max_element(work_.begin(), work_.end());
        auto max_load = std::accumulate(work_.begin(), work_.end(), 0);
        
        while(min_load != max_load)
        {
            auto guess_load  = min_load + (max_load - min_load)/2;
            auto current_load = 0;
            int guess_workers = 1;
            for(auto work : work_)
            {
                if(current_load + work > guess_load )
                {
                    ++guess_workers;
                    current_load = work;
                }
                else
                {
                    current_load += work;
                }
            }
            if(guess_workers > workers_ )
            {
                min_load = guess_load + 1;
            }
            else
            {
                max_load = guess_load;
            }
        }
        return min_load;
    }
    vector<vector<int>> GetLoads()
    {
        vector<vector<int>> result;
        int min_load = MinimizeLoad();
        vector<int> this_part;
        auto current_load = 0;
        for(auto work :work_)
        {
            if(work + current_load > min_load)
            {
                current_load = work;
                result.push_back(move(this_part));
                this_part.push_back(work);
            }
            else
            {
                current_load += work;
                this_part.push_back(work);
            }
        }
        if(this_part.size())
        {
            result.push_back(move(this_part));
        }
        return result;
    }
private:
    vector<int> work_;
    int workers_;
};
class Manachers
{
public:
    Manachers(string input):input_(input)
    {
        
    }
    string GetLongestPalindrome()
    {
        //i - c = c - imirr -> i_mirr = 2*c -i
        vector<size_t> dp(input_.size());
        size_t c = 0;
        size_t r = 0;
        for(size_t i = 1; i != input_.size(); ++i)
        {
            
            if(i < r)
            {
                size_t i_mirror = c*2 - i;
                dp[i] = min(r-i,dp[i_mirror]);
            }
            size_t r_max = i + dp[i] + 1;
            if(r_max == input_.size())
            {
                continue;
            }
            size_t l_max = ((dp[i] + 1) > i == false) ? i - (dp[i] + 1) : input_.size();
            while(input_.c_str()[r_max] == input_.c_str()[l_max])
            {
                ++dp[i];
                if(r_max > r)
                {
                    c = i;
                    r = r_max;
                }
                r_max =  i + dp[i] + 1;
                l_max = ((dp[i] + 1) > i == false) ? i - (dp[i] + 1) : input_.size();
                if(r_max == input_.size() || l_max == input_.size())
                {
                    continue;
                }
            }
            
        }
        return string{input_.begin() + c - (r-c) ,input_.begin() + r + 1};
    }
private:
    string input_;
};
int MaxSumDiscontArray(vector<int> input)
{
    vector<int> dp(input.size() + 1);
    dp[1] = input[0];
    for(int i = 2; i != input.size() + 1; ++i)
    {
        dp[i] = max(dp[i-1],input[i-1] + dp[i-2]);
    }
    dp.back();
    return 0;
}
/*
 x
 */
class KthElemIn2Arrays
{
public:
    KthElemIn2Arrays(vector<int> left, vector<int> right, size_t k):left_{left},right_{right},k_{k}
    {
        
    }
    int Compute()
    {
        
        auto left_size = left_.size();
        auto right_size = right_.size();
        if(right_size < left_size)
        {
            swap(left_,right_);
            swap(left_size,right_size);
        }
        auto end = min(left_size,k_);
        auto begin = 0;
        while(begin != end)
        {
            auto mid = begin + (end - begin)/2;
            
        }
        return 0;
    }
    
private:
    vector<int> left_,right_;
    size_t k_;
};
//


/*
 city 0 - 2 slots
 city 1=  2 slots
 1     {city0,city1}
 2     {city0,city1}
 3     {city0,city1}
 4     {city0,city1}
 
     0
 1
 
 tuple {user0,user1,city_slot0,cityslot1}
 user0,user1, 0,2
 user0,user1,1,1
 
 user3
 user3,city0,city1
 */
class AssignTravelers
{
public:
    AssignTravelers(vector<pair<int,int>> costs):costs_{move(costs)}
    {
        
    }
    pair<vector<size_t>,vector<size_t>> Compute()
    {
        pair<vector<size_t>,vector<size_t>> result;
        return result;
    }
private:
    
    
    vector<pair<int,int>> costs_;
};

class Dijsktra
{
public:
    struct Edge
       {
           size_t from_;
           size_t to_;
           double weight_;
       };
    
    using EdgeList = vector<Edge>;

    Dijsktra(const map<size_t,EdgeList>& graph,size_t src, size_t dest, int hops = numeric_limits<int>::max())
    :graph_{graph},
    predecessor_(graph.size()),
    distance_to_(graph.size(),numeric_limits<double>::max()),
    hop_to_(graph.size(),numeric_limits<int>::max()),
    min_pq_{EdgeDistanceComparator{*this}},
    max_hops_{hops + 1}
    {
        src_ = src;
        dest_ = dest;
    }
    
    void DoDijktra()
    {
        distance_to_[src_] = 0.0;
        hop_to_[src_] = 0;
        predecessor_[src_] = src_;
        min_pq_.insert(src_);
        while(!min_pq_.empty())
        {
            auto min_vertex = *min_pq_.begin();
            min_pq_.erase(min_pq_.begin());
            Relax(min_vertex);
        }
    }
    vector<size_t> GetPath()
    {
        vector<size_t> result;
        if(distance_to_[dest_] == numeric_limits<double>::max())
        return result;
        result.push_back(dest_);
        for(auto pred = predecessor_[dest_];pred != dest_;pred = predecessor_[dest_])
        {
            result.push_back(pred);
            dest_ = pred;
        }
        return result;
    }
private:
    struct EdgeDistanceComparator
    {
        EdgeDistanceComparator(Dijsktra& dijkstra):dijkstra_{dijkstra}
        {
            
        }
        bool operator()(size_t left, size_t right)
        {
            if(left == right)
                return false;
            if(dijkstra_.distance_to_[left] == dijkstra_.distance_to_[right])
                return left < right;
            return dijkstra_.distance_to_[left] < dijkstra_.distance_to_[right] ;
        }
    private:
        Dijsktra& dijkstra_;
    };
    void Relax(size_t src)
    {
        for(const auto& edge : graph_[src])
        {
            if(edge.to_ == dest_)
            {
                if(distance_to_[src_] + edge.weight_ < distance_to_[edge.to_] && hop_to_[src_] < max_hops_)
                {
                    distance_to_[edge.to_] = distance_to_[src] + edge.weight_;
                    hop_to_[edge.to_] = hop_to_[src] + 1;
                    predecessor_[edge.to_] = src;
                }
            }
            if(distance_to_[src] + edge.weight_ < distance_to_[edge.to_] && hop_to_[src] + 1 < max_hops_)
            {
                distance_to_[edge.to_] = distance_to_[src] + edge.weight_;
                hop_to_[edge.to_] = hop_to_[src] + 1;
                predecessor_[edge.to_] = src;
                min_pq_.erase(edge.to_);
                min_pq_.insert(edge.to_);
            }
        }
    }
    vector<size_t> predecessor_;
    vector<double> distance_to_;
    vector<int> hop_to_;
    map<size_t,EdgeList> graph_;
    set<size_t,EdgeDistanceComparator> min_pq_;
    int max_hops_;
    size_t src_;
    size_t dest_;

};
class CircularArraySum
{
public:
    CircularArraySum(vector<int> src):src_{move(src)}
    {
        
    }
    int Compute()
    {
        int result = 0;
        ComputeLeft();
        ComputeRight();
        for(int i = 0; i != src_.size(); ++i)
        {
            result = max(result,left_[i] + right_[i]);
        }
        return max(result,RegularMax());
    }
private:
    //1,2,3,4,-5,6
    void ComputeLeft()
    {
        int running = 0;
        for(auto i : src_)
        {
            left_.push_back(max(running,running + i));
            running += i;
        }
    }
    
    void ComputeRight()
    {
        right_.resize(src_.size());
        auto i = src_.size() - 1;
        int running = 0;
        while(i)
        {
            --i;
            right_[i] = max(right_[i+1],running += src_[i+1]);
        }
    }
    int RegularMax() const
    {
        int max_result = 0;
        int running_sum = 0;
        for(auto i : src_)
        {
            running_sum += i;
            max_result = max(max_result,running_sum < 0 ? (running_sum = 0) : running_sum);
        }
        return max_result;
    }
    vector<int> src_;
    vector<int> left_;
    vector<int> right_;
};


class FlightPath
{
public:
    FlightPath(vector<vector<string>>& ticket):tickets_{ticket}
    {
        
    }
    void ConstructGraph()
    {
        map<string,int> vertices_;
        int vertex_id = 0;
        for(auto ticket : tickets_)
        {
          if(!vertices_.count(ticket[0]))
          {
              vertices_[ticket[0]] = vertex_id++;
          }
          if(!vertices_.count(ticket[1]))
          {
              vertices_[ticket[1]] = vertex_id++;
          }
        }
        src_id_ = vertices_["JFK"];
        for(auto ticket : tickets_)
        {
            auto src_id = vertices_[ticket[0]];
            auto dst_id = vertices_[ticket[1]];
            codes_[src_id] = ticket[0];
            codes_[dst_id] = ticket[1];
            graph_[src_id].insert(FlightEdge{dst_id,ticket[1]});
        }
    }
    vector<string> BuildIternary()
    {
        ConstructGraph();
        DoEulerianPath(src_id_);
        reverse(iternary_.begin(), iternary_.end());
        return iternary_;
    }
private:
    void DoEulerianPath(int vertex)
    {
        while(!graph_[vertex].empty())
        {
            auto flight = graph_[vertex].begin();
            auto dest = flight->dest_;
            graph_[vertex].erase(flight);
            string text = codes_[vertex] + " to " + codes_[dest];
            cout << text << endl;
            DoEulerianPath(dest);
        }
        iternary_.push_back(codes_[vertex]);
    }
    struct FlightEdge
    {
        int dest_;
        string code_;
        bool operator < (const FlightEdge& edge) const
        {
            return code_.compare(edge.code_) < 0 ;
        }
    };
    map<int,string> codes_;
    vector<vector<string>>& tickets_;
    map<int,multiset<FlightEdge>> graph_;
    int src_id_;
    vector<string> iternary_;
    vector<int> iternary_idx_;
};
#include <sstream>
class UniquePaths
{
public:
    UniquePaths(int m , int n):m_{m}, n_{n},
    dp_(m,(vector<vector<string>>(n)))
    {
        
    }
    vector<string> ComputeAllPaths()
    {
        ComputeAllPaths(m_-1,n_-1);
        return dp_[m_-1][n_-1];
    }
private:
    void ComputeAllPaths(int i, int j )
    {
        if(i < 0 || j < 0)
            return;
        
        if(!dp_[i][j].empty())
            return;
        if(i == 0 && j ==0)
        {
            dp_[i][j].push_back(MakePath(i, j));
            return;
        }
        
        
        ComputeAllPaths(i-1,j);
        ComputeAllPaths(i,j-1);
        auto this_elem = MakePath(i, j);
        
        if(i)
        {
            for(auto path :dp_[i-1][j])
            {
                dp_[i][j].push_back(path + this_elem);
            }
        }
        if(j)
        {
            for(auto path :dp_[i][j-1])
            {
                dp_[i][j].push_back(path + this_elem);
            }
        }
        
    }
    string MakePath(int i, int j)
    {
        stringstream str;
        str << "{" << i << "," << j << "}";
        return str.str();
    }
    int m_;
    int n_;
    vector<string> paths_;
    vector<vector<vector<string>>> dp_;
};
class DecimalPowerset
{
public:
    DecimalPowerset(int n)
    {
        string result;
        DecimalHelper(n,result);
    }
    void DecimalHelper(int n, string& result)
    {
        if(n ==0)
        {
            cout << result << endl;
            return;
        }
        for(int i = 0; i != 10; ++i)
        {
            char ch = i + '0';
            result.push_back(ch);
            DecimalHelper(n-1,result);
            result.pop_back();
        }
    }
};
class Perms
{
public:
    Perms(string in):
    in_{move(in)}
    {
        Compute(in_.begin(),in_.end());
    }
    void  Compute(string::iterator cursor,string::iterator end)
    {
        if(cursor == end)
        {
            cout << in_ << endl;
            return;
        }
        Compute(cursor + 1, end);
        for(auto next = cursor + 1;next != end ; ++next)
        {
            if(*next != *cursor)
            {
                iter_swap(cursor, next);
                Compute(cursor + 1, end);
                iter_swap(cursor, next);
            }
        }
    }
private:
    string in_;
};
class KthSentence
{
public:
    KthSentence(string first_word, int k):k_{k}
    {
        GetKthSentence(first_word,1);
    }
    string GetKthSentence()const
    {
        return result_;
    }
private:
    void GetKthSentence(string prefix, float probability)
    {
        if(!k_)
            return;
        
        auto suggestions = GetSuggestions(prefix);
        
        if(suggestions.empty())
        {
            if(--k_== 0)
                result_ = prefix;
            
            return;
        }
        
        priority_queue<Suggestion,vector<Suggestion>> suggestions_queue;
        for(const auto& suggestion : suggestions)
        {
            auto next_suffix = prefix + " " + suggestion.word_;
            auto next_probability = probability * suggestion.probability_;
            suggestions_queue.push(Suggestion{next_suffix,next_probability});
        }
        
        while (!suggestions_queue.empty() && k_) {
            auto next_suggestion = suggestions_queue.top();
            GetKthSentence(next_suggestion.word_, next_suggestion.probability_);
            suggestions_queue.pop();
        }
    }
    struct Suggestion
    {
        string word_;
        float probability_;
        bool operator < (const Suggestion& suggestion) const
        {
            return probability_ < suggestion.probability_;
        }
    };

    vector<KthSentence::Suggestion> GetSuggestions(string prefix)
    {
        if(prefix == "I")
        {
            return {Suggestion{"will",0.8},Suggestion{"was",0.05}};
        }
        else if (prefix == "I will")
        {
            return {Suggestion{"run",0.8},Suggestion{"swim",0.6}};
        }
        else if (prefix == "I was")
        {
            return {Suggestion{"running",0.8},Suggestion{"swimming",0.6}};
        }
        vector<Suggestion> end_of_suggestion;
        return end_of_suggestion;
    }
    int k_;
    string result_;
};
/*
 cost[i][j] = 0 , if i == j
 for k - [i,j)
 cost is min over k = [i,j)
 cost[i][j] = cost[i][k] + cost[k+1][j] + dims_[i]*dims_[k+1]* dims_[j+1]
 
 */
class MatrixMultiplier
{
public:
    
    MatrixMultiplier(vector<pair<int, int>> matrices):costs_(matrices.size(),vector<int>(matrices.size()))
    {
        dims_.push_back(matrices[0].first);
        for(const auto& matrix: matrices)
        {
            dims_.push_back(matrix.second);
        }
        
        len_ = matrices.size();
        
    }
    
public:
    //A1,A2,A3
    int ComputeCost()
    {
        for(int width = 1 ; width != len_ ; ++width)
            for(int i = 0, j=i+width ; j < len_ ; ++i,j=i+width)
                for(int k = i ; k !=j ; ++k)
            {
                auto new_cost = costs_[i][k] + costs_[k+1][j] + dims_[i]*dims_[k+1]*dims_[j+1];
                if(costs_[i][j] == 0)
                {
                    costs_[i][j] = new_cost;
                }
                else
                {
                    costs_[i][j] = min(costs_[i][j] ,new_cost);
                }
                
            }
        return costs_[0][len_-1];
    }
    vector<int> dims_;
    vector<vector<int>> costs_;
    size_t len_;
};
/*
 12345
     3
 1 2   4 5
 */
//class TreeNode
class PhoneDial
{
public:
    PhoneDial() = default;
    vector<string> GetResult(const vector<int>& input)
    {
        DoDFS(input.cbegin(),input.cend());
        return result_;
    }
private:
    void DoDFS(const vector<int>::const_iterator cursor, const vector<int>::const_iterator end)
    {
        if(cursor == end)
        {
            result_.push_back(scratch_);
            return;
        }
        
        for(const auto& ch : key_map[*cursor])
        {
            scratch_.push_back(ch);
            DoDFS(cursor +1, end);
            scratch_.pop_back();
        }
    }
    vector<string> result_;
    string scratch_;
    static map<int,string> key_map ;
};
 map<int,string> PhoneDial:: key_map =
{
    
    {1,"ABC"},
    {2,"DEF"},
    {3,"GHI"},
    {4,"JKL"},
    {5,"MNO"},
    {6,"PQR"},
    {7,"STU"},
    {8,"VWX"},
    {9,"YZ"},
    
};
template <typename T>
class Slist;
template<typename T>
class SLink
{
    friend class Slist<T>;
public:
    SLink() = default;
    SLink(T data):data_(data)
    {}
    T Data()
    {
        return data_;
    }
    SLink& operator++()
    {
        if(next_)
        {
            *this = *next_;
            return *this;
        }
        return *this;
    }
    
    shared_ptr<SLink<T>> next_ = nullptr;
private:
    T data_;
    
};

class Dlink
{
public:
    explicit Dlink(int data):data_{data}
    {
        
    }
    shared_ptr<Dlink> prev_ = nullptr, next_ = nullptr;
    int data_;
};
class Dlist
{
public:
    Dlist() = default;
    Dlist& operator +=(int data)
    {
        tail->next_ = make_shared<Dlink>(data);
        tail->next_->prev_ =  tail;
        tail = tail->next_;
        ++size_;
        return *this;
    }
    void InOrder(const shared_ptr<Dlink>& root)
    {
        if(root == nullptr)
            return;
        InOrder(root->prev_);
        cout << root->data_ << endl;
        InOrder(root->next_);
    }
    Dlist& operator+= (const shared_ptr<Dlink>& link)
    {
        tail->next_ = link;
        tail = link;
        ++size_;
        return *this;
    }
    // A -> B -> C -> D
    // A <- B <- C <- D
    //swap A and C
    shared_ptr<Dlink> nth_elem(int n)
    {
        auto head = dummy->next_;
        while(n && head)
        {
            --n;
            head =head->next_;
        }
        return head;
    }
    friend ostream& operator << (ostream& os, const Dlist& list)
    {
        auto head = list.dummy->next_;
        os << "{ ";
        while(head)
        {
            os << head->data_ ;
            head = head->next_;
            if(head)
            {
                os << " , ";
            }
        }
        os << " }";
        return os;
    }
    shared_ptr<Dlink> ToBST()
    {
        return ToBST(0,size_,&dummy->next_);
    }
    shared_ptr<Dlink> ToBST(size_t begin, size_t end, shared_ptr<Dlink>* cursor)
    {
        if(begin == end)
        {
            return nullptr;
        }
        auto mid = begin + (end-begin)/2;
        auto left = ToBST(begin, mid,cursor);
        auto root = *cursor;
        root->prev_ = left;
        *cursor = (*cursor)->next_;
        root->next_ = ToBST(mid + 1, end,cursor);
        return root;
    }
    void Swap(const shared_ptr<Dlink>& left, const shared_ptr<Dlink>& right)
    {
        
        swap(left->prev_,right->prev_);
        swap(left->next_, right->next_);
        if(right->prev_)
        {
            right->prev_->next_ = right;
        }
        if(right->next_)
        {
            right->next_->prev_ = right;
        }
        if(left->next_)
        {
            left->next_->prev_ = left;
        }
        if(left->prev_)
        {
            left->prev_->next_ = left;
        }

    }
    void Reverse()
    {
        auto first = dummy->next_;
        auto last = tail;
        for(;;)
        {
            if(last == first || first->next_ == last)
                    return;
                
            Swap(first,last);
            swap(first,last);
            first = first->next_;
            last = last->prev_;
        }
    }
private:
    shared_ptr<Dlink> dummy = make_shared<Dlink>(0);
    shared_ptr<Dlink> tail = dummy;
    size_t size_=0;
};

template <typename T>
class Slist
{
public:
    Slist() = default;
    Slist(const shared_ptr<SLink<T>>& head)
    {
        dummy_head->next_ = head;
        tail = head;
    }
    shared_ptr<SLink<T>> get_head()
    {
        return dummy_head->next_;
    }
    Slist& operator +=(const T & data)
    {
        tail->next_ = make_shared<SLink<T>>(data);
        tail = tail->next_;
        return *this;
    }
    Slist& operator +=(shared_ptr<SLink<T>> insert)
    {
        tail->next_ = insert;
        tail = tail->next_;
        return *this;
    }
    //l1->l2->l3->l4
    void Sort()
    {
        //Sort(0,size_, )
    }
    void Splice(shared_ptr<SLink<T>> first, shared_ptr<SLink<T>> second)
    {
        tail->next_ = first;
        tail = second;
    }
    pair<shared_ptr<SLink<T>>,shared_ptr<SLink<T>>> Merge(pair<shared_ptr<SLink<T>>,shared_ptr<SLink<T>>> left, pair<shared_ptr<SLink<T>>,shared_ptr<SLink<T>>> right)
    {
        Slist sorted;
        auto left_begin = left.first;
        auto left_end = left.second->next_;
        
        auto right_begin = right.first;
        auto right_end = right.second->next_;
        
        while(left_begin != left_end)
        {
            if(right_begin == right_end)
            {
                sorted.Splice(left_begin, left.second);
            }
            if(right_begin->data_ < left_begin->data_)
            {
                sorted += right_begin;
                right_begin = right_begin->next_;
            }
            else
            {
                sorted += left_begin;
                left_begin = left_begin->next_;
            }
        }
        if(right_begin != right_end)
        {
            sorted.Splice(right_begin, right.second);
        }
        return {sorted.get_head(),sorted.tail};
    }
    pair<shared_ptr<SLink<T>>,shared_ptr<SLink<T>>> Sort(size_t begin, size_t end,shared_ptr<SLink<T>>* cursor )
    {
        if(begin == end)
            return  nullptr;
        if(end - begin == 1)
        {
            return {*cursor,*cursor};
        }
        auto mid = begin + (end - begin)/2;
        auto left = Sort(begin, mid,cursor);
        auto right = Sort(mid, end, cursor);
        return Merge(left, right);
    }
    void Reverse()
    {
        shared_ptr<SLink<T>> result = nullptr;
        auto head = dummy_head->next_;
        while(head)
        {
            auto next_part = head->next_;
            head->next_ = result;
            result = head;
            head = next_part;
        }
        dummy_head->next_ = result;
    }
    friend ostream& operator << (ostream& os,  Slist<T>& slist)
    {
        auto list_impl = slist.dummy_head->next_;
        while(list_impl)
        {
            os << list_impl->data_ << ",";
            list_impl = list_impl->next_;
        }
        return os;
    }
    shared_ptr<SLink<T>> Mid()
    {
       if(! dummy_head->next_ || !dummy_head->next_->next_)
           return dummy_head->next_;
        
        auto slow = dummy_head;
        auto fast = dummy_head;
        while(fast->next_ && fast->next_->next_)
        {
            fast = fast->next_->next_;
            slow = slow->next_;
        }
        return slow;
    }
    bool IsPalindrome()
    {
        auto mid_ = Mid();
        auto next_half = mid_->next_;
        mid_->next_ = nullptr;
        //mid_
        Slist next_list{next_half};
        next_list.Reverse();
        auto left_list = dummy_head->next_;
        auto right_list = next_list.get_head();
        while(left_list && right_list)
        {
            if(left_list->data_ != right_list->data_)
                return false;
            left_list = left_list->next_;
            right_list = right_list->next_;
        }
        return true;
    }
    void ReverseSubList(int begin, int end)
    {
        
    }
    int Length()const
    {
        auto head= dummy_head->next;
        int len = 0;
        while(head)
        {
            ++len;
            head = head->next_;
        }
        return len;
    }
    shared_ptr<SLink<T>> Advance(int k)
    {
        auto head= dummy_head->next;
        while(head && k)
        {
            head = head->next_;
            --k;
        }
        //assert k == 0
        return head;
    }
private:
    shared_ptr<SLink<T>> dummy_head = make_shared<SLink<T>>(T{});
    shared_ptr<SLink<T>> tail = dummy_head;
};
class KillManxOld
{
public:
    KillManxOld(vector<string> damageChart, vector<int> bossHealth):damageChart_{move(damageChart)},
    bossHealth_{move(bossHealth)}
    {
        max_ = bossHealth.size();
        bossMaskComplete_ = (1 << max_) - 1;
    }
    int GetCost()
    {
        int result = 0;
        priority_queue<Shot,vector<Shot>,greater<Shot>> lined_up;
        lined_up.push(Shot{0,0});
        short running_mask = 0;
        while(bossMask_ != bossMaskComplete_)
        {
            auto shot = lined_up.top();
            lined_up.pop();
                
            for(int i = 0; i != max_;++i)
            {
                    if((running_mask & (1 << i)) == 0)
                    {
                        
                    }
            }
        }
        return result;
    }

private:
    struct Shot
    {
        int cost;
        int boss_index;
        bool operator > (const Shot& that)const
        {
            return cost > that.cost;
        }
    };
    //int
        
    short bossMask_ = 0;
    short bossMaskComplete_;
    vector<string>damageChart_;
    vector<int> bossHealth_;
    size_t max_;
    
};
bool visited[(1 << 15)];

struct node {
int weapons = 0;
int shots = 0;
    bool operator > (const node& other)const
    {
        return shots > other.shots;
    }
// Define a comparator that puts nodes with less shots on top appropriate to your language
};
int leastShots(vector<string> damageChart, vector<int> bossHealth) {
priority_queue<node,vector<node>,greater<>> pq;

    pq.push(node{});

while (pq.empty() == false) {
node top = pq.top();
pq.pop();

// Make sure we don’t visit the same configuration twice
if (visited[top.weapons]) continue;
visited[top.weapons] = true;
    auto numWeapons = bossHealth.size();
// A quick trick to check if we have all the weapons, meaning we defeated all the bosses.
// We use the fact that (2^numWeapons - 1) will have all the numWeapons bits set to 1.
if (top.weapons == (1 << numWeapons) - 1)
return top.shots;

for (int i = 0; i < damageChart.size(); i++) {
// Check if we’ve already visited this boss, then don’t bother trying him again
if ((top.weapons >> i) & 1) continue;

// Now figure out what the best amount of time that we can destroy this boss is, given the weapons we have.
// We initialize this value to the boss’s health, as that is our default (with our KiloBuster).
int best = bossHealth[i];
for (int j = 0; j < damageChart.size(); j++) {
if (i == j) continue;
    const char* damaged_elem = damageChart[j].c_str();
if (((top.weapons >> j) & 1) && damaged_elem[i] != '0') {
// We have this weapon, so try using it to defeat this boss
    int shotsNeeded = bossHealth[i] / (damageChart[j][i] - '0');
    if (bossHealth[i] % (damaged_elem[i] - '0') != 0) shotsNeeded++;
    best = min(best, shotsNeeded);
}
}

// Add the new node to be searched, showing that we defeated boss i, and we used ‘best’ shots to defeat him.
    pq.push(node{top.weapons | (1 << i), top.shots + best});
}
}
    return 0;
}
/*
        0
   |    |    |
{150,  150, 150}
   |    |     |
{"070","500","140"}

*/

typedef long long ll;
class KiloManX
{
  public:
  
  size_t n;
  vector<string> damageChart;
  vector<int> bossHealth;
  map<int,ll> best;
  
  ll getBest(int done)
  {
    if (done == (1<<n)-1) return 0;
    if (best.count(done)) return best[done];
    
    ll bestCost = 12345678912345;
    for (int nb = 0; nb < n; nb++)
    if (! (done & (1 << nb)) )
    {
      int nextDone = done + (1 << nb);
      int md = 1;
      for (int j = 0; j < n; j++)
      {
          if (done & (1 << j))
              md = max(md,damageChart[j][nb]-'0');
      }
      ll cost = bossHealth[nb];
      cost += (md-1);
      cost /= md;
      auto next_cost = getBest(nextDone);
      bestCost = min (bestCost, cost + next_cost);
    }
    
    best[done] = bestCost;
    return bestCost;
  }
  
  int leastShots(vector <string> _damageChart, vector <int> _bossHealth)
  {
    damageChart = _damageChart;
    bossHealth = _bossHealth;
    n = damageChart.size();
    
    return getBest(0);
  }
};
struct GraphVertex
{
    int id_;
    int distance_;
    bool operator < (const GraphVertex& that)const
    {
        if(id_ == that.id_)
            return false;
        if(distance_ < that.distance_)
            return true;
        if(distance_ == that.distance_)
        {
            return id_ < that.id_;
        }
        return false;
    }
};
struct GraphVertexComp
{
    static map<int,int> distance;
    bool operator () (int id_left, int id_right) const
    {
        if(id_left == id_right)
            return false;
        
        if(distance[id_left] == distance[id_right])
            return id_left < id_right;
        
        return distance[id_left] < distance[id_right];
    }
};
map<int,int> GraphVertexComp::distance = {{2,5},{0,5}, {1,6}};

class InflectionPoint {
public:
    int search(const vector<int>& nums, int target) {
        if(nums.empty())
            return -1;
    
        auto inflection = inflection_point(nums);
        
        if(inflection == nums.size())
        {
            return bin_search(nums,0,inflection,target);
        }
        if(target == nums[inflection])
            return inflection;

        if(target < nums[inflection])
            return -1;
        
        if(nums[inflection -1 ] < target)
            return -1;
        
        if(target > nums[inflection])
        {
            auto res_0 =  bin_search(nums, inflection + 1, nums.size(),target);
            if(res_0 != -1)
                return res_0;
        }
        return bin_search(nums,0, inflection, target);
        
    }
    
private:
    int bin_search(const vector<int>& nums, int begin, int end, int target)
    {
        auto u_b_result = u_b(nums,begin,end,target);
        auto l_b_result = l_b(nums,begin,end,target);
        return l_b_result != u_b_result ? l_b_result : -1;
    }
    int u_b(const vector<int>& nums, int begin, int end, int target)
    {
        while(begin != end)
        {
            auto mid = begin + (end - begin)/2;
            if(target < nums[mid] )
            {
                end = mid;
                
            }
            else
            {
                begin = mid + 1;
            }
        }
        return begin;
    }
    
    int l_b(const vector<int>& nums, int begin, int end, int target)
    {
        while(begin != end)
        {
            auto mid = begin + (end - begin)/2;
            if(nums[mid] < target)
            {
                begin = mid + 1;
            }
            else
            {
                end = mid;
            }
        }
        return begin;
    }
    
    int inflection_point(const vector<int>& nums)
    {
        auto begin = 0;
        auto end = nums.size();
        while(begin != end)
        {
            auto mid = begin + (end - begin)/2;
            if(nums[mid] < nums[0])
            {
                end = mid;
            }
            else
            {
                begin = mid + 1;
            }
        }
        return begin;
    }
    
};

  struct TreeNodeX {
     int val;
      TreeNodeX *left;
      TreeNodeX *right;
      TreeNodeX() : val(0), left(nullptr), right(nullptr) {}
      TreeNodeX(int x) : val(x), left(nullptr), right(nullptr) {}
      TreeNodeX(int x, TreeNodeX *left, TreeNodeX *right) : val(x), left(left), right(right) {}
  };


class TSolution {
public:
    using MaxSum_RunningSUM =  pair<int,int>;
    int maxPathSum(TreeNodeX* root) {
        if(root == nullptr)
            return 0;
        if(root->left == nullptr && root->right == nullptr)
            return root->val;
        int current_max = std::numeric_limits<int>::min();
        maxPathSum(root,current_max);
        return current_max;
    }
private:
    int maxPathSum(TreeNodeX* root, int &current_max)
    {
        if(!root)
        {
            return 0;
        }
        auto left = max(0,maxPathSum(root->left,current_max));
        auto right = max(0,maxPathSum(root->right,current_max));
        auto this_max = max(current_max,root->val + left + right);
        current_max = max(this_max,current_max);
        return root->val + max(left,right);
    }
       
};
#include <unordered_map>
class LRUCache {
public:
    LRUCache(int capacity):capacity_{capacity} {
        
    }
    
    int get(int key) {
        
        return get_helper(key);
    }
    
    void put(int key, int value) {
        put_helper(key, value);
    }
private:
    int get_helper(int key)
    {
        if(cache_.count(key) == 0)
           return -1;
        lru_list_.erase(cache_[key].lru_list_entry_);
        lru_list_.push_front(key);
        cache_[key].lru_list_entry_ = lru_list_.cbegin();
        return cache_[key].val_;
    }
    void put_helper(int key, int value)
    {
        if(cache_.count(key))
        {
            lru_list_.erase(cache_[key].lru_list_entry_);
        }
        else if(cache_.size() == capacity_)
        {
            EvictLRUKey();
        }

        lru_list_.push_front(key);
        cache_[key] = CacheElem{value,lru_list_.cbegin()};
    }
    void EvictLRUKey()
    {
        cache_.erase(lru_list_.back());
        lru_list_.pop_back();
    }
    struct CacheElem
    {
        int val_;
        list<int>::const_iterator lru_list_entry_;
    };
    unordered_map<int,CacheElem> cache_;
    list<int> lru_list_;
    int capacity_;
};
class SolutionWeight {
public:
    int lastStoneWeight(vector<int>& stones) {
        for(auto stone : stones)
        {
            stones_.push(stone);
        }
        return stoneit();
    }
    string largestTimeFromDigits(vector<int>& A) {
        
        //1,2,3,4
        sort(A.begin(),A.end());
        int i = 4;
        char result[6] {0};
        
        while(i)
        {
            --i;
            if(A[i] < 4)
            {
                if(i > 0)
                {
                    if(A[i] < 3)
                    {
                        result[0] = A[i] + '0';
                        result[1] = A[i-1]+ '0';
                    }
                    else
                    {
                        result[0] = A[i-1] + '0';
                        result[i] = A[i] + '0';
                    }
                    break;
                }
            }
        }
        switch (i) {
            case 3:
            {
                result[3] = A[1] +  '0';
                result[4] = A[0] + '0';
                result[2] = ':';
                break;
            }
            case 2:
            {
                if(A[3] < 6)
                {
                    result[3] = A[3] + '0';
                    result[4] = A[0] + '0';
                    result[2] = ':';
                }
                else if(A[3] < 10)
                {
                    result[3] = A[0] + '0';
                    result[4] = A[3] + '0';
                    result[2] = ':';
                }
                break;
            }
            case 1:
            {
                if(A[3] < 6)
                {
                    result[3] = A[3] + '0';
                    result[4] = A[2] + '0';
                    result[2] = ':';
                }
                else if(A[3] < 10 && A[2] < 6)
                {
                    result[3] = A[2] + '0';
                    result[4] = A[3] + '0';
                    result[2] = ':';
                }
            }
            default:
                break;
        }
        if(result[2] == ':')
            return result;
        
        return "";
        
    }
private:
    int stoneit()
    {
        while(!stones_.empty())
        {
            auto stone_0 = stones_.top();
            stones_.pop();
            
            if(stones_.empty())
                return stone_0;
            auto stone_1 = stones_.top();
            stones_.pop();
            if(stone_1 != stone_0)
            {
                stones_.push(stone_0 - stone_1);
            }
                
        }
        return 0;
    }
    priority_queue<int> stones_;
};
namespace LongestPar
    {
    class Solution {
    public:
        
        int longestValidParentheses(string s)
        {
            stack<char> tracker;
            int last_end = -1;
            int last_max_ = 0;
            
            int current_start = -1;
            int cursor = -1;
            int running_max = 0;
            for(auto ch: s)
            {
                ++cursor;
                if(ch == ')')
                {
                    if(tracker.empty())
                        continue;
                    if(tracker.top() == '(')
                    {
                        tracker.pop();
                        running_max +=2;
                        if(tracker.empty())
                        {
                            if(last_end != -1)
                            {
                              if(last_end + 1 == current_start)
                              {
                                  last_end = cursor;
                                  last_max_ += running_max;
                              }
                            }
                            else
                            {
                                last_end = cursor;
                                last_max_ = running_max;
                            }
                            max_ = max(max_,last_max_);
                        }
                        
                    }
                    else
                    {
                        running_max = 0;
                        stack<char> empty{};
                        tracker.swap(empty);
                    }
                }
                else if(ch == '(')
                {
                    if(tracker.empty())
                    {
                        current_start = cursor;
                    }
                    tracker.push(ch);
                }
            }
            max_ = max(max_,running_max);
            return max_;
        }
    private:
        int max_ = 0;
    };
    }
class RollingWindow
{
    //"barfoofoobarthefoobarman"
    //aa,aa
    //aaaaaa
public:
    RollingWindow(const string& src, const string& text)
    {
        result_idx_ = Search(src,text);
        
    }
    vector<int> findAnagrams(string s, string p)
    {
        Search(s,p);
        return result;
    }
    
    int Search(const string& src, const string& text)
    {
        auto src_len = src.size();
        auto text_len = text.size();
        
        if(src_len < text_len)
            return src_len;
        
        auto remaining = text.size();
        
        for(auto ch : text)
        {
            ++sliding_window_[ch].first;
        }
        
        for(int i = 0; i != text_len; ++i)
        {
            if(sliding_window_.count( src[i]))
            {
                auto& [expected, found] = sliding_window_[src[i]];
                ++found;
                if(found > expected == false)
                {
                    --remaining;
                    if(!remaining)
                    {
                        result.push_back(i + 1 - text_len);
                    }
                }
            }
        }
        
        for(int j = text_len; j != src_len; ++j )
        {
            if(sliding_window_.count(src[j - text_len]))
            {
                auto& [expected, found] = sliding_window_[src[j-text_len]];
                --found;
                if(found < expected)
                {
                    ++remaining;
                }
            }
            if(sliding_window_.count(src[j]))
            {
                auto& [expected, found] = sliding_window_[src[j]];
                ++found;
                if(found > expected == false)
                {
                    --remaining;
                    if(!remaining)
                    {
                        result.push_back(j + 1 - text_len);
                    }
                }
            }
        }
        return src_len;
    }
private:
    vector<int> result;
    
    using Match = pair<int,int>;
    unordered_map<char, Match> sliding_window_;
    RollingWindow(const vector<string>& texts):texts_(texts)
    {
        text_len_ = texts_[0].size();
        for(int i = 0; i!= text_len_; ++i)
        {
            power_ = i ? (power_ * kBase_) % kMod_ : 1;
        }
        for(const auto& text:texts_)
        {
           ++histogram_[text];
        }
    }
    int ComputeHash(const string& text)
    {
        int result = 0;
        for( auto ch : text)
        {
            result = ((result * power_) + ch) % kMod_;
        }
        return result;
    }
    int power_ = 1;
    const int kBase_ = 26;
    const int kMod_ = 997;
    int text_len_=0;
    map<string, int> histogram_;
    map<string, int> window_;
    const vector<string> texts_;
    int result_idx_=-1;
};
namespace ThreeSum
{
    class Solution {
    public:
        int threeSumClosest(vector<int>& nums, int target) {
            int result = numeric_limits<int>::max();
            auto len = nums.size();
            sort(nums.begin(),nums.end());
            
            for( int j = 0; j != len; ++j)
            {
                result = min(result,threeSumClosesthelper(nums,target - nums[j],j));
            }
            
            return result;
        }
        private:
        int threeSumClosesthelper(vector<int>& nums, int target, int skip_index)
        {
            auto len = nums.size();
            int result = numeric_limits<int>::max();
            size_t begin = 0, end = len -1;
            while(begin != end )
            {
                if(begin == skip_index || end == skip_index)
                    continue;
                
                result = min(result,abs(nums[begin] + nums[end] - target));
                if(nums[begin] + nums[end] > target)
                {
                    --end;
                }
                else if(nums[begin] + nums[end] < target)
                {
                    ++begin;
                }
                else
                    return target;
            }
            return result;
        }
    };
}
/*
 0
 3
8 25
 {3,8,25}
 {0,3,3}
 */
class SvcDag
{
public:
    SvcDag(vector<int> svc_parents, vector<int> svc_children)
    {
        auto max_size = svc_parents.size();
        for(int idx = 0; idx != max_size; ++idx)
        {
            graph_[svc_parents[idx]].push_back(svc_children[idx]);
            ++indegree_[svc_children[idx]];
        }
        
    }
    

    vector<int> TopologicalSort(int src)
    {
        topo_queue_.push(src);
        while(!topo_queue_.empty())
        {
            auto top = topo_queue_.front();
            topo_queue_.pop();
            for(auto child :graph_[top] )
            {
                if(--indegree_[child] == 0)
                {
                    result_.push_back(child);
                    topo_queue_.push(child);
                }
            }
        }
        return result_;
    }
    void TopSortDFS(int src)
    {
        visited_[src] = true;
        for(auto child : graph_[src])
        {
            if(!visited_[child])
            {
                TopSortDFS(child);
            }
        }
        result_.push_back(src);
        
    }
    vector<int> result_;
    queue<int> topo_queue_;
    unordered_map<int, vector<int>> graph_;
    unordered_map<int, int> indegree_;
    unordered_map<int, bool> visited_;
};
class Heap
{
public:
    Heap(vector<int> input):input_{input}
    {
        size_ = input_.size();
        Heapify();
    }
    void Heapify()
    {
        int size = size_;

        size /= 2;
        
        while( size )
        {
            --size;
            sink(size);
        }
    }
    int pop()
    {
        if(size_)
        {
            --size_;
            swap(input_[0],input_[size_]);
            sink(0);
            return input_[size_];
        }
        throw "";
    }
    bool empty()const
    {
        return size_ == 0;
    }
private:
    void sink(int index)
    {
        int min_index = -1;
        while((min_index  =  (index*2 -1)) < size_)
        {
            if(min_index + 1 < size_)
            {
                min_index = input_[min_index + 1] < input_[min_index ] ? min_index + 1 : min_index;
            }
            if(input_[index] < input_[min_index ] == false)
            {
                swap(input_[index], input_[min_index]);
                index = min_index;
            }
            else
            {
                break;
            }
        }
        
    }
    int size_;
    vector<int> input_;
};

template <class RandIter>
void InsertSort(RandIter first, RandIter last)
{
  if(first == last || first == --last)
    return;
   for(auto cursor = first + 1; cursor != last + 1; ++cursor)
        for(auto insert_pos = cursor; insert_pos != first && *insert_pos < *(insert_pos -1);--insert_pos)
            iter_swap(insert_pos, insert_pos - 1);
        
}
                
template <class RandIter>
bool NextPerm(RandIter first, RandIter last)
{
        if(first == last || first == --last)
            return false;
        
        auto end = last;
        while(end != first && *(end -1) < *end == false)
        {
            --end;
        }
        if(end == first)
        {
            reverse(end,last + 1);
            return false;
        }
        auto greater = last;
        while(*(end -1) < *greater == false)
        {
            --greater;
        }
        iter_swap(end -1, greater);
        reverse(end,last + 1);
        
        return true;
}
struct TNode
{
    int data_;
    TNode* left_, * right_;
};

TNode* InvertTree(TNode* root)
{
    if(!root)
        return nullptr;
    if(!root->left_)
    {
        return root;
    }
    auto result = InvertTree(root->left_);
    root->left_->right_ = root;
    root->left_->left_ = root->right_;
    root->left_ = root->right_ = nullptr;
    return result;
}
namespace Tarjan
{
    class Solution {
    public:
        vector<vector<int>> criticalConnections(int n, vector<vector<int>>& connections)
        {
            BuildGraph(n,connections);
            discovery_time_.resize(n);
            low_link_.resize(n);
            DFS_Tarjan(0,0);
            return result_;
        }
    private:
        int DFS_Tarjan(int node, int parent)
        {
            low_link_[node] = discovery_time_[node] = ++discovery_timer_;
            for(auto adj_node : graph_[node])
            {
                if(adj_node == parent)
                {
                    continue;
                }
                //Zero discovery time means unvisited
                if( discovery_time_[adj_node] ==  0 )
                {
                    auto low_link = DFS_Tarjan(adj_node,node);
                    if(discovery_time_[node] < low_link )
                    {
                        result_.push_back({node,adj_node});
                    }
                    else
                    {
                        low_link_[node] = min(low_link_[node],low_link);
                    }
                }
                else
                {
                    low_link_[node] = min(low_link_[node],discovery_time_[adj_node]);
                }
            }
            return low_link_[node];
        }
        void BuildGraph(int nodes, vector<vector<int>>& connections)
        {
            graph_.resize(nodes);
            for(const auto& edge: connections)
            {
                graph_[edge[0]].push_back(edge[1]);
                graph_[edge[1]].push_back(edge[0]);
            }
        }
        vector<vector<int>> graph_;
        vector<int> discovery_time_;
        vector<int> low_link_;
        int discovery_timer_ = 0;
        vector<vector<int>> result_;
    };
}
namespace WordLadder
    {
    class Solution {
    public:
        vector<vector<string>> findLadders(string beginWord, string endWord,  vector<string>& wordList) {
            
            int idx = 0;
            for(const auto& word : wordList)
            {
                word_idx_[word] = idx++;
            }
            if(word_idx_.count(beginWord) ==0)
            {
                wordList.push_back(beginWord);
                word_idx_[beginWord] = idx;
            }
            
            BuildGraph(wordList);
            if(word_idx_.count(endWord) ==0)
                return vector<vector<string>>{};
            beginWord_ = beginWord;
             
            auto start_list = GetNeighbors(beginWord,wordList);
            auto dest = word_idx_[endWord];
            return DoBFS(word_idx_[beginWord],dest,wordList);
        }
    private:
        vector<vector<string>> DoBFS(int src, int dest, const vector<string>& wordList)
        {
            vector<vector<string>> result_;
            vector<int> dist_to(word_idx_.size(),-1);
            
            vector<int> bfs_tree(word_idx_.size(),-1);
            queue<int> bfs_queue;
            bfs_queue.push(src);
            dist_to[src] = 0;
            while(!bfs_queue.empty())
            {
                auto parent = bfs_queue.front();
                bfs_queue.pop();
                
                for(auto adj : graph_[parent])
                {
                    if(dist_to[adj] != -1)
                    {
                        if(adj == dest && dist_to[parent] + 1 == dist_to[dest])
                        {
                            auto path =  BuildPath(bfs_tree, parent,src,dest);
                            result_.push_back(BuildResult(src,path, wordList));
                        }
                        continue;
                    }
                        
                    bfs_tree[adj] = parent;
                    if(adj == dest)
                    {
                        auto path =  BuildPath(bfs_tree, parent,src,dest);
                        result_.push_back(BuildResult(src,path, wordList));
                    }
                    dist_to[adj] = dist_to[parent] + 1;
                    bfs_queue.push(adj);
                }
                
            }
            
            return result_;
        }
        vector<int> BuildPath(vector<int>& bfs_tree, int leaf,int src,int dest)
        {
            vector<int> path;
            path.push_back(dest);
            path.push_back(leaf);
            int parent = -1;
            while(( parent = bfs_tree[leaf]) != -1)
            {
                path.push_back(parent);
                leaf = parent;
            }
            //path.push_back(src);
            reverse(path.begin(),path.end());
            return path;
        }
        vector<string> BuildResult(int src, const vector<int>& path,const vector<string>& wordList)
        {
            
            vector<string> word_ladder;
            
                         for(auto elem: path)
                         {
                             word_ladder.push_back(wordList[elem]);
                         }
                         //result.push_back(move(word_ladder));
            return word_ladder;
        }
        void BuildGraph(const vector<string>& wordList)
        {
            //["hot","dot","dog","lot","log","cog"]
            auto word_count = wordList.size();
            graph_.resize(word_count);
            for(auto word_idx = 0; word_idx != word_count; ++ word_idx)
                for(auto adj_word= word_idx +1; adj_word != word_count; ++adj_word)
                {
                    if(IsNeighbor(wordList[word_idx], wordList[adj_word]))
                    {
                        graph_[word_idx].push_back(adj_word);
                        graph_[adj_word].push_back(word_idx);
                    }
                }
            
        }
        
        bool IsNeighbor(const string& first, const string& second)
        {
            auto size = first.size();
            bool mismatchFound = false;
            for(auto idx =0; idx != size; ++idx)
            {
                if(first.at(idx) != second.at(idx))
                {
                    if(mismatchFound)
                        return false;
                    mismatchFound = true;
                }
            }
            return mismatchFound;
        }
        vector<int> GetNeighbors(const string& input,const vector<string>& wordList)
        {
            vector<int> result;
            for(const auto& word: wordList)
            {
                if(IsNeighbor(word, input))
                {
                    result.push_back(word_idx_[word]);
                }
            }
            return result;
        }
        string beginWord_;
        vector<vector<int>> graph_;
        unordered_map<string, int> word_idx_;
        
    };
    }
class SplitTree
{
    struct tree_node
    {
        tree_node(int data):data_(data)
        {
            
        }
        int data_;
        int sum_ = 0;
        tree_node* left_ = nullptr,*right_ = nullptr;
        
    };

public:
    SplitTree(tree_node* root):root_(root)
    {
        sum_ = BuildSum(root);
    }
public:
    void SplitIt(tree_node* root)
    {
        if(root == nullptr)
            return ;
        //break left link
        auto right_sum =  root->data_ + root->right_->sum_;
        auto left_split_sum = sum_ - right_sum;
        
        auto left_split_product = left_split_sum * right_sum;
        max_sum_product_ = max(max_sum_product_ , left_split_product);
        
        //break right link
        auto left_sum = root->data_ + root->left_->sum_;
        auto right_split_sum = sum_ - left_sum;
        auto right_split_prod = left_sum * right_split_sum;
        max_sum_product_ = max(max_sum_product_ , right_split_prod);
        SplitIt(root->left_);
        SplitIt(root->right_);
    }
    int BuildSum(tree_node* root)
    {
       if(root == nullptr)
           return 0;
        root->sum_ = BuildSum(root->left_);
        root->sum_ += BuildSum(root->right_);
        root->sum_ += root->data_;
        return root->sum_;
    }
    tree_node* root_;
    int sum_ = 0;
    int max_sum_product_ = numeric_limits<int>::min();
    
};
namespace SubsetSumSolution
    {
    class SubsetSum
    {
    public:
        SubsetSum(const vector<int>& input):input_(input)
        {
            //sort(input_.begin(),input_.end());
        }
        vector<int> GetSubsetX(int target)
        {
            
            vector<vector<int>> dp(target + 1,vector<int>(input_.size(),-1));
            dp[input_[0]][0] = 0;
            for(int i = 1; i != target +1; ++i)
            {
                for(int elem_idx = 1; elem_idx != input_.size(); ++ elem_idx)
                {
                    if(i < input_[elem_idx] == false)
                    {
                        if(dp[i - input_[elem_idx]][elem_idx -1] != -1 || (i == input_[elem_idx]))
                        {
                            dp[i][elem_idx] = elem_idx;
                        }
                        else
                        {
                            dp[i][elem_idx] = dp[i][elem_idx-1];
                        }
                    }
                }
            }
            vector<int> result;
            if(dp[target][input_.size()-1] != -1)
            {
                auto i = input_.size()-1;
                while(target)
                {
                    if(dp[target][i] != dp[target][i-1])
                    {
                        result.push_back(input_[i]);
                        target -= input_[i];
                    }
                    --i;
                }
            }
            return result;
        }
        vector<int> GetSubset(int target)
        {
            vector<int> result;
            vector<int> dp(target + 1,-1);
            
            auto input_size = input_.size();
            for(int i = 1; i !=target + 1; ++i)
            {
                for(int elem = 0; elem != input_size; ++elem)
                {
                    if(input_[elem] > i == false)
                    {
                        if((dp[i - input_[elem]] != -1) || (i - input_[elem] == 0))
                        {
                            dp[i] = elem;
                            if(i == target)
                            {
                                result = BuildResult(target, dp);
                            }
                        }
                    }
                }
            }
            result = BuildResult(target, dp);
            return result;
            
        }
    private:
        vector<int> BuildResult(int target, const vector<int>& dp) const
        {
            vector<int> result;
            while(target && dp[target] != -1)
            {
                result.push_back(input_[dp[target]]);
                target -= input_[dp[target]];
            }
            return result;
        }
         vector<int> input_;
    };
    }
namespace DyanamicProg
{
/*
 1 2 3
 4 5 6
 7 8 9
 * 0 #
 */
    class EnumPhoneNumbers
    {
    public:
        EnumPhoneNumbers(int num_digits, unordered_map<int, vector<int>> dialPad):num_digits_{num_digits},
        dialPad_{move(dialPad)}
        {
            
        }
        int GetCount()
        {
            vector<vector<int>> dp(2,vector<int>(10));
            for(auto& elem:dp[0])
            {
                elem = 1;
            }
            for(int i = 1; i != num_digits_; ++i)
            {
                
            }
            int result = 0;
            return result;
        }
    private:
        int num_digits_;
        unordered_map<int, vector<int>> dialPad_;
    };
}
int main(int argc, const char * argv[])
{
    SubsetSumSolution::SubsetSum sssum{{2,5,3,7}};
    auto ss_res = sssum.GetSubset(10);
    
    
    WordLadder::Solution sln_ladder;
    
    /*
     "red"
     "tax"
     ["ted","tex","red","tax","tad","den","rex","pee"]
     */
    //{"a","b","c"};
    vector<string> word_list{"ted","tex","red","tax","tad","den","rex","pee"};//{"hot","dot","dog","lot","log","cog"};
    sln_ladder.findLadders("red", "tax", word_list);
    //sln_ladder.findLadders("a", "c", {"a","b","c"});
    
   //[[0,1],[1,2],[2,0],[1,3]]
    vector<vector<int>> tarjan_edges = {{0,1}, {1,2},{2,0},{1,3}};
    Tarjan::Solution tarjan;
    auto crit_sections = tarjan.criticalConnections(4, tarjan_edges);
                
    vector<int> permin{5,4,3,2,10,5};
    InsertSort(permin.begin(), permin.end());
    copy(permin.begin(), permin.end(),ostream_iterator<int>{cout,","});
    while(NextPerm(permin.begin(), permin.end()))
    {
        copy(permin.begin(), permin.end(),ostream_iterator<int>{cout,","});
        cout << endl;
    }
    cout << "real perm : " << endl ;
    
    while(next_permutation(permin.begin(), permin.end()))
    {
        copy(permin.begin(), permin.end(),ostream_iterator<int>{cout,","});
        cout << endl;
    }
     /*
      0
     1 2
    3 4
      */
    vector<int> jumper = {2,3,1,1,4};
    JumpGame jg{jumper};
    auto min_jumps = jg.BFS();

    vector<int> svc_ch{3,8,25};
    vector<int>svc_p {0,3,3};
    
    
    SvcDag svcdag{svc_p,svc_ch};
    svcdag.TopSortDFS(0);
    
    Heap hp{svc_ch};
    vector<int> sorted;
    while(!hp.empty())
    {
        sorted.push_back(hp.pop());
    }
    
    vector<int> myvec{1,2,3};
    
    swap(myvec[0],myvec[2]);
    
    RollingWindow rw{"cbaebabacd", "abc"};
    
    TreeNodeX tn{-2};
    TreeNodeX tn_l{-1};
    tn.left = &tn_l;
    
    TSolution sln;
    sln.maxPathSum(&tn);
    
    InflectionPoint pt;
    auto found = pt.search({3,1}, 3);
    
    CircularArraySum cas{{1,2,3,-5,10}};
    auto circular_max = cas.Compute();
    
    set<int,GraphVertexComp> map_comp;
    map_comp.insert(0);
    map_comp.insert(1);
    map_comp.insert(2);
    
    map_comp.erase(1);
    
    set<GraphVertex> nodes;
    nodes.insert(GraphVertex{0,5});
    nodes.insert(GraphVertex{1,5});
    
    nodes.erase(GraphVertex{1,5});
    
    //nodes.insert
    vector<string> damage_chart{"070","500","140"};
    vector<int> boss_health{150,150,150};
    KiloManX km;
    auto least_shots = km.leastShots(damage_chart, boss_health);
    least_shots = leastShots(damage_chart,boss_health);
    Dlist dlist;
    
    for(auto elem :{1,2,3,4,5})
    {
        dlist += elem;
    }
    auto bst = dlist.ToBST();
    auto left = dlist.nth_elem(0);
    auto right = dlist.nth_elem(3);
    dlist.Reverse();
    //dlist.Swap(left, right);
    cout << dlist;
    
    Slist<int> my_list;
    for(auto elem :{1,2,1})
    {
        my_list += elem;
    }
    auto is_palindrome = my_list.IsPalindrome();
    //Slist rev{ my_list.Mid()};
    //cout << rev;
    
    my_list.Reverse();
    cout << my_list;
    Slist<int> even_list, odd_list;
    auto head = my_list.get_head();
    while(head)
    {
        even_list+= head->Data();
        head = head->next_;
        if(head)
        {
            odd_list += head->Data();
           head = head->next_;
        }
    }
    cout << my_list << endl << "even:" << even_list << endl << "odd:" << odd_list;
    StockMaster ms{{1,4,2},2};
    auto max_profit = ms.MaxProfitWithCooldown();
    
    //[2,3,1,1,4]
    
    vector<pair<int, int>> mtxes = {{30,35},{35,15},{15,5},{5,10},{10,20},{20,25}};
    MatrixMultiplier mtx{move(mtxes)};
    auto cost = mtx.ComputeCost();
    
    Paranthesize ps{3};
    KthSentence kth{"I",5};
    kth.GetKthSentence();
    Perms pr{"abc"};
    
    UniquePaths paths{4,5};
    auto path_result = paths.ComputeAllPaths();
    //[["JFK","KUL"],["JFK","NRT"],["NRT","JFK"]]
    //{ {"JFK","KUL"}, {"JFK", "NRT"}, {"NRT","JFK"}}
    vector<vector<string>> sched= {{"EZE","AXA"},{"TIA","ANU"},{"ANU","JFK"},{"JFK","ANU"},{"ANU","EZE"},{"TIA","ANU"},{"AXA","TIA"},{"TIA","JFK"},{"ANU","TIA"},{"JFK","TIA"}};
    
    //{ {"JFK","SFO"}, {"JFK", "ATL"}, {"SFO","ATL"}, {"ATL","JFK"}, {"ATL","SFO"} };
    FlightPath fp{sched};
    auto iter = fp.BuildIternary();
    using EdgeList = vector<Dijsktra::Edge>;
    map<size_t,EdgeList> flights
    {
        {0,{ {0,1,100.0},{0,2,500.0} } },
        {1,{{1,2,200.0}}},
    };
    Dijsktra dj{flights,0,2,1};
    dj.DoDijktra();
    auto path = dj.GetPath();
    
    auto disc_array = MaxSumDiscontArray({2,5,3,1,7});
    Manachers mch{"abacabacabb"};
    auto pali_max = mch.GetLongestPalindrome();
    WorkBalancer wb{vector<int>{ 10, 20, 60, 50, 30, 40},3};
    auto parts = wb.GetLoads();
    
    vector<int> width = {9236, 7065, 2283, 506, 6432, 9812, 3133, 1397, 7052, 3729,
        2556, 9954, 1367, 6440, 5141, 3091, 2879, 1346, 7080, 1036,
        7503, 7775, 433, 7579, 6520, 2287, 1971, 3879, 1725, 8200,
        1830, 2774, 3850, 7509, 8531, 7493, 1511, 9399, 9679, 2124,
        791, 3432};
    vector<int> speed = {956, 799, 481, 194, 993, 444, 571, 986, 815, 910,
        98, 847, 650, 487, 419, 434, 410, 812, 374, 751,
        307, 134, 134, 955, 758, 73, 932, 360, 135, 588,
        218, 936, 674, 494, 157, 556, 881, 292, 851, 890,
        886, 912};
    
    OptimumTravel opt{width,speed,9756};
    
    auto best_time = opt.ComputeTime();
    AnaGramm ag{"search text", "ttxe"};
    auto ag_pos = ag.Find();
    makeAnagramEx("cde","abc");
    jumpingOnClouds(vector<int>{0,1,0,0,0,1,0});
    
    CountSubsequences cs{"rabbbit","rabbit"};
    auto count_result = cs.Count();
    ConcurrentVector<int> cv;
    cv.push_back(5);
    //cv[0] = 0;
    AFunc(cv[0]);
    
    EqualSumBackTrack equalSumsBT{vector<int>{1,2,4,2,3,6,1,5},6};
    equalSumsBT.DoBackTrack(0);
    KnapSack ks
    {
        vector<pair<int,int>>{
        {2,4},
        {3,2},
        {4,10}
        }
    };
    ks.Compute(5);
    RegExp regExp{"ab*"};
    regExp.Match("abc");
    
    SplitFairly spf{{10,11,5,6,10}};
    spf.Compute();
    PalindromicDecompositon pdecmp{"sliril"};
    
    
   vector<vector<int>> grid =
    {
        {0,0,0},
        {0,1,0},
        {0,0,0}
    }
   /* {{0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},{0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,0,0},{1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,1},{0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0},{0,0,0,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0},{1,0,1,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0},{0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,0,0,0,0,0},{0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},{1,0,1,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,0,1,0,0,0,1,0,1,0,0,0,0,1},{0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,0,0,0,0},{0,1,0,1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,0,0,0,0,0},{0,1,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,1},{1,0,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0},{0,0,1,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,1,1,0,1,0,0,0,0,1,1},{0,1,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,1,0,1},{1,1,1,0,1,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0,0},{0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0,1,1},{0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,1,0,0,0}}
    */
    ;
    Solution soln;
    auto max_paths = soln.uniquePathsWithObstacles(grid);
    
    
    lcs("rabbbit","rabbit");
    auto distinct = CountDistinct("rabbit","rabbbit");
    MinWindow minw{"ADOBECODEBANC", "ABC"};
    minw.GetMin();
    ms.GetTransactions();
    LCS lcs_{"cable", "dable"};
    auto max_lcs = lcs_.Compute();
    
    vector<int> zigs =  { 1,3,2,4,1,3,1,4,5,2,2,1,4,2,2};
    WaterMark wm{zigs};
    wm.CalculateArea();
    
    LongPal pl{"lirils"};
    auto long_pal = pl.LongPalinDrome();
    auto hex = HextToString(255);
    auto res1 = AnagramPalindrome("lirill");
    
   /*
    Person('P1', 2001, 2005),
    Person('P2', 2005, 2009),
    Person('P3', 2003, 2020),
    Person('P3', 2002, 2002),
    */
    string ip = "ababcbacadefegdehijhklij";
    MaxPartition mp{ip};
    mp.PartitionUnique();
    vector<int> people = {3,5,3,4};
    
    auto boats = MaxBoats(people,5);
    
    vector<vector<int>> grid_gate =
    {
        {std::numeric_limits<int>::max(),-1,0,std::numeric_limits<int>::max()},
        {std::numeric_limits<int>::max(),std::numeric_limits<int>::max(),std::numeric_limits<int>::max(),-1},
        {std::numeric_limits<int>::max(),-1,std::numeric_limits<int>::max(),-1},
        {0,-1,std::numeric_limits<int>::max(),std::numeric_limits<int>::max()}
    };
    GridLocked gl{grid_gate};
    gl.FillGrid();
    string vow = "hello";
    ReverseVowels(vow);
    
    auto can_pali = CanPalindrome("lirill");
    vector<string> text = {"apple",
        {"mango"},
        {"apple"},
        {"banana"},
        {"Berries"}
    };
    vector<string> kw = {"apple",
        {"banana"}
    };
    FindInWindow fw(text,kw);
    fw.GetWindow();
    vector<TrafficVolume> traffics = {{0,1.3},{1,0.0},{2,2.5},{3,3.7},{4,0.0}, {5,1.4},{6,2.6},
        {7,0.0},
        {8,2.2},
        {9,1.7},
        {10,0},
        {11,0.0},
        {12,0.0},
        {13,0.0},
        {14,1.7}
        
    };
    SlidingWindowMax slw{ traffics};
    auto max_sizes = slw.GetMax(3);
    vector<int> input_array = {5,4,2,1,3,3,3};
    auto res = QuickSelect(input_array.begin(), input_array.end(),6);
    
    vector<Event> registry = {
        Event{2001,Event::Type::Arrived},
        Event{2005,Event::Type::Arrived},
        Event{2003,Event::Type::Arrived},
        Event{2002,Event::Type::Arrived},
        Event{2005,Event::Type::Departed},
        Event{2009,Event::Type::Departed},
        Event{2020,Event::Type::Departed},
        Event{2002,Event::Type::Departed},
        
    };
    NPR npr(registry);
    auto max_alive = npr.MaxPopulation();
    
        auto infection = FindInflection({2,3,4,5,6,-1,0,1});
    
    pivotIndex({1,7,3,6,5,6});
/*    1. even
 4 /2 = 2
 n/2 , n/2 -1
 0,1,2,3
 
 
    //2,3,5, 6,8,9,10,11,12
    //1,4,7,13  ,14,15
 
     
    //1,2,3,    4,5,6,
    //7,8,9,10, 11, 12,13,14,15
    //
 
 */
    auto array_0 = {1, 3, 6, 30};
    ContiguousSum cont_sum{array_0};
    cont_sum.CanSum(31);
    
    HistoGram hg{zigs};
    auto area = hg.GetResult();
    {
        string input = "lracecartt";
        LongestPalindromic pl(input);
        auto result = pl.GetResult();
    }
    std::vector<TimeElement> time_series =
    {
        {2,2},
        {3,3},
        {1,3},
        {1,2},
        {1,1},
        {1,0},
    };
    Chronicler chr{4};
    chr.Normalize(time_series);
    
;
    
    auto max_zig = ZigZag(zigs);
    

    auto sqrt = BinSq(300);
    
    int lc_max = lcs("harry", "sally");
    auto count = stepPerms(5);
    EnumSubstrings enumerator{"kkkk"};
    enumerator.Enum();
    int anagrams = enumerator.CountAnagrams();
    
   vector<vector<int>> grid_vec=
    {
        {1, 1, 0, 0},
        {0, 1, 1, 0},
        {0, 0, 1, 0},
        {1, 0, 0, 0}
    };
    Regional reg{grid_vec};
    auto max_region = reg.GetMax();
    
    SuffixArray arr{L"banana"};
    MinCoins mc{{1,5,6,9}};
    int targetSum = 9;
    auto mincoins = mc.DoMinFast(targetSum);
    mc.DoMinCoin(targetSum);
    auto result = mc.GetMinCoins(targetSum);
    auto ways = mc.DoWays(11);
    
    //auto result = mc.GetMinCoins(targetSum);
   // std::ostream_iterator<int> os{cout,","};
    //std::copy(result.cbegin(),result.cend(),os);
    //auto ways = mc.DoWays(11);
    
    string src = "eskimos";
    string target = "smokes";
    EditDistance ed{src,target};
    
    // insert code here...
    std::vector<std::vector<size_t>> sudoku =
    {
        {0,0,0, 3,0,0, 0,9,0},
        {0,0,0, 0,1,9, 0,6,0},
        {0,0,0, 0,5,0, 3,0,8},
        
        {0,0,5, 0,0,3, 8,0,0},
        {4,9,0, 0,7,0, 0,1,5},
        {0,0,6, 1,0,0, 9,0,0},
        
        {9,0,2, 0,6,0, 0,0,0},
        {0,5,0, 4,8,0, 0,0,0},
        {0,6,0, 0,0,7, 0,0,0}
    };
    Sudoku sdk{sudoku};

    return 0;
}
//5,2,3,4,6,1

//1,2,3,4,5,6
//0,6,2,3,4,1,5
//i = 1, j = 5

/*
 1 - 6 -> i = 1,
 6 - 5 -> i =
 1,6,5,1
 */

/*
 Welcome to Facebook!

 This is just a simple shared plaintext pad, with no execution capabilities.

 When you know what language you would like to use for your interview,
 simply choose it from the dropdown in the top bar.

 Enjoy your interview!


 Question:  list, k. return whether we can find contigious sequence of integer form list that sum to k

 example: list = [1, 3, 4, 30], k=7, return true. (3, 4)
          list = [1, 3, 6, 30], k=7, return false.
  #include <vector>
   #include <queue>
   using std::vector;
 using std::queue;

 queue: (1, 3, 4)

 list = [1, 1, 1, 1, 1, 30, 3], k=31
   
 list = [1, 3], k = 4
   
 bool CanSum(vector<int> input, int k)
 {
     queue<int> window;
     int current_sum = 0;
   
     for(auto elem:input)
     {
       if(current_sum == k)
         return true;
       
       while(current_sum > k && !queue.empty())
       {
         
         current_sum -= window.front();
         window.pop();
       }
       else
       {
         current_sum += elem;
         window.push(elem);
         
         
       }
       
     }
     
     return current_sum == k;
 }


 Q2: list , k. return kth smallest element from list. time compleixty O(nlogk) or better
 example: list = [3, 5, 1, 4, 2] k = 2, return 2


 int GetKth(vector<int> input, int k)
 {
   if(input.size() < k)
     throw std::exception("Invlaid input");
  
   std::priority_queue<int> pq;
   for(int i = 0; i != k; ++i)
   {
     pq.push(input[i]);
   }
   
   for(int j = k; j != input.size())
   {
     if(input[j]) > pq.top)
       continue;
     pq.push(input[j]);
     pq.pop();
   }
   
   return pq.top();
 }
 

 **/
/*
 #include <bits/stdc++.h>

 using namespace std;



 /*
  
 #include <bits/stdc++.h>
 #include <limits>

 using namespace std;



 /*
  * Complete the 'minimumWindowSubstring' function below.
  *
  * The function is expected to return a STRING.
  * The function accepts following parameters:
  *  1. STRING fullString
  *  2. STRING chars
  */

 

