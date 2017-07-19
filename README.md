<h4>Computation of periodic orbits for non-autonomous systems</h4>

These routines compute the initial conditions for a periodic orbit of a non-autonomous with periodic dependence on time:
<P></P>
<DIV ALIGN="CENTER" CLASS="mathdisplay"><!-- MATH
 \begin{equation*}
\dot{x}=f(x,t),\quad f(x,t+T)=f(x,t).
\end{equation*}
 -->
<TABLE CLASS="equation*" CELLPADDING="0" WIDTH="100%" ALIGN="CENTER">
<TR VALIGN="MIDDLE">
<TD NOWRAP ALIGN="CENTER"><SPAN CLASS="MATH"><IMG
 WIDTH="246" HEIGHT="32" ALIGN="MIDDLE" BORDER="0"
 SRC="img1.png"
 ALT="$\displaystyle \dot{x}=f(x,t),\quad f(x,t+T)=f(x,t).$"></SPAN></TD>
<TD NOWRAP CLASS="eqno" WIDTH="10" ALIGN="RIGHT">
&nbsp;&nbsp;&nbsp;</TD></TR>
</TABLE></DIV>
<BR CLEAR="ALL"><P></P>
That is, the goal is to compute a point <SPAN CLASS="MATH"><IMG
 WIDTH="21" HEIGHT="14" ALIGN="BOTTOM" BORDER="0"
 SRC="img2.png"
 ALT="$ x^*$"></SPAN> such  that, for a given
<SPAN CLASS="MATH"><IMG
 WIDTH="17" HEIGHT="30" ALIGN="MIDDLE" BORDER="0"
 SRC="img3.png"
 ALT="$ t_0$"></SPAN>,
<P></P>
<DIV ALIGN="CENTER" CLASS="mathdisplay"><!-- MATH
 \begin{equation*}
\varphi(t_0+qT;x)=x.
\end{equation*}
 -->
<TABLE CLASS="equation*" CELLPADDING="0" WIDTH="100%" ALIGN="CENTER">
<TR VALIGN="MIDDLE">
<TD NOWRAP ALIGN="CENTER"><SPAN CLASS="MATH"><IMG
 WIDTH="130" HEIGHT="32" ALIGN="MIDDLE" BORDER="0"
 SRC="img4.png"
 ALT="$\displaystyle \varphi(t_0+qT;x)=x.$"></SPAN></TD>
<TD NOWRAP CLASS="eqno" WIDTH="10" ALIGN="RIGHT">
&nbsp;&nbsp;&nbsp;</TD></TR>
</TABLE></DIV>
<BR CLEAR="ALL"><P></P>
The method used consists on performing a Newton method to solve the equation
<P></P>
<DIV ALIGN="CENTER" CLASS="mathdisplay"><!-- MATH
 \begin{equation*}
S_{t_0}(x)-x=0,
\end{equation*}
 -->
<TABLE CLASS="equation*" CELLPADDING="0" WIDTH="100%" ALIGN="CENTER">
<TR VALIGN="MIDDLE">
<TD NOWRAP ALIGN="CENTER"><SPAN CLASS="MATH"><IMG
 WIDTH="110" HEIGHT="32" ALIGN="MIDDLE" BORDER="0"
 SRC="img5.png"
 ALT="$\displaystyle S_{t_0}(x)-x=0,$"></SPAN></TD>
<TD NOWRAP CLASS="eqno" WIDTH="10" ALIGN="RIGHT">
&nbsp;&nbsp;&nbsp;</TD></TR>
</TABLE></DIV>
<BR CLEAR="ALL"><P></P>
where <SPAN CLASS="MATH"><IMG
 WIDTH="26" HEIGHT="30" ALIGN="MIDDLE" BORDER="0"
 SRC="img6.png"
 ALT="$ S_{t_0}$"></SPAN> is the stroboscopic map (time <SPAN CLASS="MATH"><IMG
 WIDTH="24" HEIGHT="30" ALIGN="MIDDLE" BORDER="0"
 SRC="img7.png"
 ALT="$ qT$"></SPAN> return map):
<P></P>
<DIV ALIGN="CENTER" CLASS="mathdisplay"><!-- MATH
 \begin{equation*}
S_{t_0}(x)=\varphi(t_0+qT;x)
\end{equation*}
 -->
<TABLE CLASS="equation*" CELLPADDING="0" WIDTH="100%" ALIGN="CENTER">
<TR VALIGN="MIDDLE">
<TD NOWRAP ALIGN="CENTER"><SPAN CLASS="MATH"><IMG
 WIDTH="159" HEIGHT="32" ALIGN="MIDDLE" BORDER="0"
 SRC="img8.png"
 ALT="$\displaystyle S_{t_0}(x)=\varphi(t_0+qT;x)$"></SPAN></TD>
<TD NOWRAP CLASS="eqno" WIDTH="10" ALIGN="RIGHT">
&nbsp;&nbsp;&nbsp;</TD></TR>
</TABLE></DIV>
<BR CLEAR="ALL"><P></P>
The program needs the gsl libraries, and can be compiled using
<DIV ALIGN="CENTER">
g++ -o main main.c rk78.c lu.c memory.c -lgsl -lgslcblas

</DIV>
