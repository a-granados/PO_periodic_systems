Computation of periodic orbits for non-autonomous systems

These routines compute the initial conditions for a periodic orbit of a non-autonomous with periodic dependence on time:
<DIV ALIGN="CENTER" CLASS="mathdisplay"><!-- MATH
 \begin{equation*}
\dot{x}=f(x,t),\quad f(x,t+T)=f(x,t)
\end{equation*}
 -->
<TABLE CLASS="equation*" CELLPADDING="0" WIDTH="100%" ALIGN="CENTER">
<TR VALIGN="MIDDLE">
<TD NOWRAP ALIGN="CENTER"><SPAN CLASS="MATH"><IMG
 WIDTH="241" HEIGHT="32" ALIGN="MIDDLE" BORDER="0"
 SRC="img1.png"
 ALT="$\displaystyle \dot{x}=f(x,t),\quad f(x,t+T)=f(x,t)$"></SPAN></TD>
<TD NOWRAP CLASS="eqno" WIDTH="10" ALIGN="RIGHT">
&nbsp;&nbsp;&nbsp;</TD></TR>
</TABLE></DIV>

The method used is to solve a Newton method to solve the equation
<DIV ALIGN="CENTER" CLASS="mathdisplay"><!-- MATH
 \begin{equation*}
S_{t_0}(x)-x=0
\end{equation*}
 -->
<TABLE CLASS="equation*" CELLPADDING="0" WIDTH="100%" ALIGN="CENTER">
<TR VALIGN="MIDDLE">
<TD NOWRAP ALIGN="CENTER"><SPAN CLASS="MATH"><IMG
 WIDTH="105" HEIGHT="32" ALIGN="MIDDLE" BORDER="0"
 SRC="img2.png"
 ALT="$\displaystyle S_{t_0}(x)-x=0$"></SPAN></TD>
<TD NOWRAP CLASS="eqno" WIDTH="10" ALIGN="RIGHT">
&nbsp;&nbsp;&nbsp;</TD></TR>
</TABLE></DIV>
where
for the stroboscopic map:
<DIV ALIGN="CENTER" CLASS="mathdisplay"><!-- MATH
 \begin{equation*}
S_{t_0}(x)=\varphi(t_0+qT;x)
\end{equation*}
 -->
<TABLE CLASS="equation*" CELLPADDING="0" WIDTH="100%" ALIGN="CENTER">
<TR VALIGN="MIDDLE">
<TD NOWRAP ALIGN="CENTER"><SPAN CLASS="MATH"><IMG
 WIDTH="152" HEIGHT="32" ALIGN="MIDDLE" BORDER="0"
 SRC="img3.png"
 ALT="$\displaystyle S_{t_0}(x)=\varphi(t_0+qT;x)$"></SPAN></TD>
<TD NOWRAP CLASS="eqno" WIDTH="10" ALIGN="RIGHT">
&nbsp;&nbsp;&nbsp;</TD></TR>
</TABLE></DIV>
is the stroboscopic map (time qt return map).
The solution to this equation provides a point such that
<DIV ALIGN="CENTER" CLASS="mathdisplay"><!-- MATH
 \begin{equation*}
\varphi(t_0+qT;x)=x
\end{equation*}
 -->
<TABLE CLASS="equation*" CELLPADDING="0" WIDTH="100%" ALIGN="CENTER">
<TR VALIGN="MIDDLE">
<TD NOWRAP ALIGN="CENTER"><SPAN CLASS="MATH"><IMG
 WIDTH="125" HEIGHT="32" ALIGN="MIDDLE" BORDER="0"
 SRC="img4.png"
 ALT="$\displaystyle \varphi(t_0+qT;x)=x$"></SPAN></TD>
<TD NOWRAP CLASS="eqno" WIDTH="10" ALIGN="RIGHT">
&nbsp;&nbsp;&nbsp;</TD></TR>
</TABLE></DIV>
and x is the initial condition for a qT periodic orbit.
