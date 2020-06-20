#!/usr/bin/env python
# coding: utf-8

# In[31]:

print('This code generates the Forman-Ricci curvature for the COVID-19 data provided by the John Hopkins University:')

print('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv'
)

print('The Forman-ricci curvature for simple graphs is also implemented and can be found at:  ')

print('https://github.com/saibalmars/GraphRicciCurvature')


print('Packages required:')
print('-pandas')
print('-numpy')
print('-networkx')
print('-matplotlib')
print('-seaborn')
print('-ssl')

print('In case of issues, please, contact us: danillo.dbs16@gmail.com')
print('Process started. The whole process might take a several time, depending on your computation power.')

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


# In[32]:


# The following Forman-Ricci curvature for networks can be found at https://github.com/saibalmars/GraphRicciCurvature


# In[33]:


def formanCurvature(G, verbose=False):
    """
     Compute Forman-ricci curvature for all nodes and edges in G.
         Node curvature is defined as the average of all it's adjacency edge.
     :param G: A connected NetworkX graph, unweighted graph only now, edge weight will be ignored.
     :param verbose: Show detailed logs.
     :return: G: A NetworkX graph with Forman-Ricci Curvature with node and edge attribute "formanCurvature"
     """

    # Edge forman curvature
    for (v1, v2) in G.edges():
        if G.is_directed():
            v1_nbr = set(list(G.predecessors(v1)) + list(G.successors(v1)))
            v2_nbr = set(list(G.predecessors(v2)) + list(G.successors(v2)))
        else:
            v1_nbr = set(G.neighbors(v1))
            v1_nbr.remove(v2)
            v2_nbr = set(G.neighbors(v2))
            v2_nbr.remove(v1)
        face = v1_nbr & v2_nbr
        # G[v1][v2]["face"]=face
        prl_nbr = (v1_nbr | v2_nbr) - face
        # G[v1][v2]["prl_nbr"]=prl_nbr

        G[v1][v2]["formanCurvature"] = len(face) + 2 - len(prl_nbr)
        if verbose:
            print("Source: %s, target: %d, Forman-Ricci curvature = %f  " % (v1, v2, G[v1][v2]["formanCurvature"]))

    # Node forman curvature
    for n in G.nodes():
        fcsum = 0  # sum of the neighbor Forman curvature
        if G.degree(n) != 0:
            for nbr in G.neighbors(n):
                if 'formanCurvature' in G[n][nbr]:
                    fcsum += G[n][nbr]['formanCurvature']

            # assign the node Forman curvature to be the average of node's adjacency edges
            G.nodes[n]['formanCurvature'] = fcsum / G.degree(n)
        if verbose:
            print("node %d, Forman Curvature = %f" % (n, G.nodes[n]['formanCurvature']))
    if verbose:
        print("Forman curvature computation done.")

    return G


# In[34]:


import ssl
ssl._create_default_https_context = ssl._create_unverified_context
JH=pd.read_csv('https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv')


# In[35]:


import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


# In[36]:


from dateutil.parser import parse


# In[37]:


import datetime


# In[38]:


countries_JH=JH['Country/Region'].unique().tolist()
dates_JH=list(JH)[4:]
N=[]
for c in countries_JH:
    L=JH[JH['Country/Region']==c][dates_JH].sum().values.tolist()
    N.append(L)
JH_df=pd.DataFrame(np.array(N).T).rename(columns={i:countries_JH[i] for i in range(len(countries_JH))})    


# In[39]:


Jh_nc=[]
for c in countries_JH:
    L=JH_df[c].to_list()

    L1=[]
    L1.append(L[0])
    for i in range(1,len(L)-1):
        L1.append([L[i+1]-L[i] if L[i+1]-L[i]>=0 else 0][0])
    Jh_nc.append(L1)    


# In[40]:


JH_nc=pd.DataFrame(np.array(Jh_nc).T).rename(columns={i:countries_JH[i] for i in range(len(countries_JH))})   


# In[41]:


Df=JH_nc.copy()


# In[42]:


def graph(t0,t1,e):
    

    data=np.array(Df[t0-1:t1].corr().fillna(0))
    M=[[data[i][j] if data[i][j]>1.0001-e else 0 for i in range(len(data))] for j in range(len(data))]

    G=nx.from_numpy_array(np.array(M))
    for n in G.nodes():
        if G.has_edge(n,n):
            G.remove_edge(n,n)
            
    return G


# In[43]:


def e(t0,t1,step=0.01):
    e=0
    G=graph(t0,t1,e)
    n=nx.number_connected_components(graph(t0,t1,2))
    while nx.number_connected_components(G)>n:
        e+=step
        G=graph(t0,t1,e)
        #print(e)
    return e    


# In[44]:


import datetime
from datetime import datetime


# In[45]:


Df=JH_nc.copy()

time=7

Forman=[]
for i in range(1,len(Df)-time+2):

    E=e(i,i+time-1)
    G=graph(i,i+time-1,E)
    formanCurvature(G)
    value=np.mean(list(nx.get_edge_attributes(G,'formanCurvature').values()))
    Forman.append(value)


# In[46]:


time=7
epi_sum=[]
for i in range(1,len(Df)-time+2):
    epi_sum.append(Df[i:i+time-1].sum().sum())


# In[47]:


dates=list(JH)[4:]


# In[48]:


dates = [datetime.strptime(i, '%m/%d/%y').strftime('%d/%b') for i in dates]


# In[49]:


moving=[dates[i]+str('-')+dates[i+time-1] for i in range(len(dates)-time)]


# In[50]:


ticks=[]
moving_ticks=[]
for i in range(len(moving)):
    if i%2!=0:
        moving_ticks.append(moving[i])
        ticks.append(i)


# In[51]:



fig=plt.figure(figsize=(30,10))
fig.text(0.1, 0.4, 'COVID worldwide vs. Forman-Ricci Curvature',fontsize=30, ha='center', va='center', rotation='vertical')
#plt.title('Time Window: 7 days')

plt.subplot(211)
plt.plot(epi_sum,lw=2,color='crimson',marker='.',label='Total new cases \n (last update: '+dates[-1]+'/2020)')
plt.xticks(ticks,['' for i in ticks])

plt.yticks(size=20)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,1))
plt.legend(framealpha=0,fontsize=20)
plt.tick_params(width=2,length=4)
plt.margins(x=0)

plt.subplot(212)
plt.plot(Forman,lw=2,color='green',marker='.',label='Mean Forman-Ricci Curvature')
plt.plot([0 for i in Forman],color='gray',ls='dashed')
plt.xticks(ticks,moving_ticks)
plt.xticks(rotation='90',size=20)
plt.xticks(size=20)
plt.yticks(size=20)
plt.tick_params(width=2,length=4)

plt.legend(framealpha=0,fontsize=20)

plt.margins(x=0)

color_index=moving_ticks.index('11/Mar-17/Mar')
plt.gca().get_xticklabels()[color_index].set_color("red")
import seaborn
seaborn.despine(top=True)
plt.savefig('Ricci_on_COVID_new_cases.png',dpi=300,bbox_inches='tight')


# In[23]:


Df=JH_df.copy()
#m=50
time=7
Forman=[]
now = datetime.now()
for i in range(1,len(Df)-time+2):
    E=e(i,i+time-1)
    G=graph(i,i+time-1,E)
    formanCurvature(G)
    value=np.mean(list(nx.get_edge_attributes(G,'formanCurvature').values()))
    Forman.append(value)


# In[24]:


epi_cum=[]
for i in range(time-1,len(JH_df)):
    epi_cum.append(JH_df[i:i+1].sum().sum())


# In[25]:


ticks=[i for i in range(len(dates[time-1:]))]


# In[26]:


shift_dates=dates[time-1:]


# In[27]:


ticks=[]
label_ticks=[]
for i in range(len(shift_dates)):
    if i%2!=0:
        ticks.append(i)
        label_ticks.append(shift_dates[i])


# In[30]:



fig=plt.figure(figsize=(30,10))
fig.text(0.09, 0.44, 'COVID worldwide vs. Forman-Ricci Curvature',fontsize=30, ha='center', va='center', rotation='vertical')


plt.subplot(211)
plt.plot(epi_cum,lw=2,color='crimson',marker='.',label='Cumulative cases \n (last update: '+dates[-1]+'/2020)')
plt.xticks(ticks,['' for i in ticks])

plt.yticks(size=20)
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,1))
plt.legend(framealpha=0,fontsize=20)
plt.tick_params(width=2,length=4)
plt.margins(x=0)

plt.subplot(212)
plt.plot(Forman,lw=2,color='green',marker='.',label='Mean Forman-Ricci Curvature')
plt.plot([0 for i in Forman],color='gray',ls='dashed')
plt.xticks(ticks,label_ticks)
plt.xticks(rotation='90',size=20)
plt.xticks(size=20)
plt.yticks(size=20)
plt.tick_params(width=2,length=4)



plt.legend(framealpha=0,fontsize=20)

plt.margins(x=0)
color_index=label_ticks.index('11/Mar')
plt.gca().get_xticklabels()[color_index].set_color("red")
import seaborn
seaborn.despine(top=True)
plt.savefig('Ricci_on_COVID_cumulative_cases.png',dpi=300,bbox_inches='tight')

print('Process termitaded. Figures saved with the following names: ')
print('Ricci_on_COVID_new_cases.png')
print('Ricci_on_COVID_cumulative_cases.png')
# In[ ]:




