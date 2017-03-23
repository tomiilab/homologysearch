class Coor:
    def __init__(self,figsize,l_r,b_t,hor_ver,row_col):
        fig_wid,fig_hei=figsize
        self.col_wid=fig_wid+sum(l_r)
        self.col_hei=fig_hei+sum(b_t)
        self.total_wid=fig_wid*row_col[1]+sum(l_r)+hor_ver[0]*(row_col[1]-1)
        self.total_hei=fig_hei*row_col[0]+sum(b_t)+hor_ver[1]*(row_col[0]-1)
        self.left=l_r[0]/self.total_wid
        self.right=1-l_r[1]/self.total_wid
        self.bottom=b_t[0]/self.total_hei
        self.top=1-b_t[1]/self.total_hei
        self.wspace,self.hspace=hor_ver
    def figsize(self):
        return {'figsize':(self.total_wid,self.total_hei)}
    def adjust(self):
        return {'wspace':self.wspace,'hspace':self.hspace,'left':self.left,'right':self.right,'top':self.top,'bottom':self.bottom}
