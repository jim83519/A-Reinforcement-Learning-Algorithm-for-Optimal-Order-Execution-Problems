
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 10 16:54:50 2018

@author: 10707305
"""

import numpy as np
import pandas as pd
import time

begin=time.time()
np.random.seed(1)

N_STATES = 3   # 1维世界的宽度
STATE=[1900,2800,3700]
remain_shares=100
ACTIONS = np.arange(0,remain_shares+1,1)     # 探索者的可用动作
EPSILON = 0.9   # 贪婪度 greedy
ALPHA = 0.01     # 学习率
GAMMA = 0.9    # 奖励递减值,就是discount rate
MAX_EPISODES = 21210000   # 最大回合数
DELTA_T=1/(N_STATES-1)
RHO=3.8
SPREAD=0.02
transition_prob=[[0.85,0.1375,0.0125],[0.1392405,0.6075949,0.2531646],[0,0.2625,0.7375]]

def build_q_table(n_states, actions):
    table = pd.DataFrame(
        np.zeros((n_states, len(actions))),     # q_table 全 0 初始
        columns=actions,    # columns 对应的是行为名称
    )
    return table

      

# 在某个 state 地点, 选择行为
def choose_action(state, q_table, remain_shares):
    state_actions = q_table[state, range(remain_shares+1)]  # 选出这个 state 的所有 action 值
    index = (state_actions==0)
    actions=np.arange(0,remain_shares+1,1)
    if (np.random.uniform() > EPSILON) or (sum(index) == len(state_actions)):  # 非贪婪 or 或者这个 state 还没有探索过
        action_name = np.random.choice(actions)
    else:
        action_name = state_actions.argmax()    # 贪婪模式
    return action_name

def get_env_feedback(S, A, D, state):
    # This is how agent will interact with the environment
    if  S == N_STATES - 1:   # terminate
        S_ = 'terminal'        
        R = -A*(D+(A/(2*state)))
        D = np.exp(-RHO*DELTA_T)*D+np.exp(-RHO*DELTA_T)*(1/state)*A
    else:
        S_ = S + 1
        R = -A*(D+(A/(2*state)))        
        D = np.exp(-RHO*DELTA_T)*D+np.exp(-RHO*DELTA_T)*(1/state)*A
    
    return S_, R, D

def rl():
    q_table_0 = np.array(build_q_table(len(STATE), ACTIONS))  # 初始 q table (*2為兩個狀態 將不同狀態的Q值都寫出來)
    q_table_1 = np.array(build_q_table(len(STATE), ACTIONS))  # 初始 q table (*2為兩個狀態 將不同狀態的Q值都寫出來)
    q_table_2 = np.array(build_q_table(len(STATE), ACTIONS))  # 初始 q table (*2為兩個狀態 將不同狀態的Q值都寫出來)
    
    for episode in range(MAX_EPISODES):     # 回合
        S_index_0=np.asscalar(np.random.choice(len(STATE),1))   #數列為[,)的型式，決定初始狀態, 可在最後設機率
        state = STATE[S_index_0]
        spread = SPREAD
        S = 0 #第一步為0開始 (0,1,2)         
        x_0 = remain_shares
        A_0 = choose_action(S_index_0, q_table_0, x_0)
        S_, R_0, spread= get_env_feedback(S, A_0, spread, state)
        x_1 = x_0-A_0
        
        if(state==STATE[0]):
            p=transition_prob[0]
        if(state==STATE[1]):
            p=transition_prob[1]
        if(state==STATE[2]):
            p=transition_prob[2]
        S_index_1=np.asscalar(np.random.choice(len(STATE),1,p=p)) #轉換狀態                
        S = S_
        state = STATE[S_index_1]
        A_1 = choose_action(S_index_1, q_table_1, x_1)
        S_, R_1, spread= get_env_feedback(S, A_1, spread, state)
        x_2 = x_1-A_1
        
        if(state==STATE[0]):
            p=transition_prob[0]
        if(state==STATE[1]):
            p=transition_prob[1]
        if(state==STATE[2]):
            p=transition_prob[2]
        S_index_2=np.asscalar(np.random.choice(len(STATE),1,p=p)) #轉換狀態
        S = S_
        state = STATE[S_index_2]
        A_2 = x_2
        S_, R_2, spread= get_env_feedback(S, A_2, spread, state)
        
        
        q_predict = q_table_0[S_index_0, A_0]    # 估算的(状态-行为)值
        q_target = R_0 + GAMMA * q_table_1[S_index_1, A_1] + GAMMA * GAMMA * q_table_2[S_index_2, A_2]
        q_table_0[S_index_0, A_0] += ALPHA * (q_target - q_predict)
        
        q_predict = q_table_1[S_index_1, A_1]    # 估算的(状态-行为)值
        q_target = R_1 + GAMMA * q_table_2[S_index_2, A_2]
        q_table_1[S_index_1, A_1] += ALPHA * (q_target - q_predict)
        
        q_predict = q_table_2[S_index_2, A_2]    # 估算的(状态-行为)值
        q_target = R_2
        q_table_2[S_index_2, A_2] += ALPHA * (q_target - q_predict)
        

        if episode % 10000 == 0:
            print(episode//10000)
        #print(episode)
        
    return q_table_0, q_table_1, q_table_2

if __name__ == "__main__":
    q_table_0, q_table_1, q_table_2 = rl()
    print('\r\nQ-table:\n')
    print(q_table_0, q_table_1, q_table_2)
    
    
status=pd.DataFrame()
for i in range(len(STATE)):
    for j in range(len(STATE)):
        for k in range(len(STATE)):
            temp = pd.DataFrame([[i,j,k]])
            status = status.append(temp)    
ans=pd.DataFrame()
for i in range(len(status)):
    s = pd.Series(status.iloc[i,:])
    D = SPREAD
    x_0 = remain_shares
    
    u_0 = q_table_0[s[0],:].argmax()
    x_0 = x_0-u_0
    
    u_1 = q_table_1[s[1],:(x_0+1)].argmax()
    x_0 = x_0-u_1
    u_2 = x_0
    temp=[u_0,u_1,u_2]
    temp1=pd.DataFrame({"first":[temp[0]],"second":[temp[1]],"third":[temp[2]]})
    ans=pd.concat([ans,temp1])
print(ans)

end=time.time()
print(end-begin)

ans.to_csv("test_dis_prob_2121_spread002.csv")