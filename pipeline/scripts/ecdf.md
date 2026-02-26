%model selection

% ---------------- BRCA ----------------
\begin{figure}[htbp]
  \centering
  \includegraphics[width=\linewidth]{Thesis_writeup/Figures/ecfd_plots/ecdf_breast_cancer.pdf}
\caption{Empirical cumulative distribution of cell line similarity ranks across
breast cancer clinical patient samples.
For each breast cancer clinical patient sample ($n = 818$), all cell lines were
ranked by tumour--cell line transcriptomic similarity (rank 1 = highest similarity).
Curves show the ECDF of the ranks achieved by the top breast cancer cell line
candidates (KPL-1, MCF-7, EFM-192A, MDA-MB-468-B), i.e.\ the fraction of patient
samples for which a given cell line attains rank $\leq x$.
The vertical dashed line denotes the top-$k$ threshold ($k = 10$), and the horizontal
dashed line marks 0.5 (median fraction).
A non-breast cancer cell line is included as a negative control (grey).
Curves that rise more steeply at low ranks indicate models that more consistently
achieve high transcriptomic similarity across patient samples.}
\label{fig:brca_ecdf_model_selection}
\end{figure}


% ---------------- NBL ----------------
\begin{figure}[htbp]
  \centering
  \includegraphics[width=\linewidth]{Thesis_writeup/Figures/ecfd_plots/ecdf_neuroblastoma.pdf}
\caption{Empirical cumulative distribution of cell line similarity ranks across
neuroblastoma clinical patient samples.
For each neuroblastoma clinical patient sample ($n = 244$), cell lines were ranked
by tumour--cell line transcriptomic similarity (rank 1 = highest similarity).
Curves show the ECDF of the ranks achieved by the top neuroblastoma cell line
candidates (NBLS, IMR32, CHP126, KELLY), i.e.\ the fraction of patient samples
for which each cell line attains rank $\leq x$.
The vertical dashed line denotes the top-$k$ threshold ($k = 10$), and the horizontal
dashed line marks 0.5 (median fraction).
A non-neuroblastoma cell line is included as a negative control (grey).
Curves that rise earlier indicate candidates that more consistently place near
the top of patient-specific rankings.}
\label{fig:nbl_ecdf_model_selection}
\end{figure}

% ---------------- RBL ----------------
\begin{figure}[htbp]
  \centering
  \includegraphics[width=\linewidth]{Thesis_writeup/Figures/ecfd_plots/ecdf_retinoblastoma.pdf}
\caption{Empirical cumulative distribution of cell line similarity ranks across
retinoblastoma clinical patient samples.
For each retinoblastoma clinical patient sample ($n = 66$), cell lines were ranked
by tumour--cell line transcriptomic similarity (rank 1 = highest similarity).
Curves show the ECDF of the ranks achieved by the top retinoblastoma cell line
candidates (RBL-14, RBL-18, RBL-20), i.e.\ the fraction of patient samples for
which a given cell line attains rank $\leq x$.
The vertical dashed line denotes the top-$k$ threshold ($k = 10$), and the horizontal
dashed line marks 0.5 (median fraction).
A non-retinoblastoma cell line is included as a negative control (grey), which
consistently ranks beyond position 10 across patient samples, in contrast to
the retinoblastoma cell lines that achieve median ranks of 2--4.}
\label{fig:rbl_ecdf_model_selection}
\end{figure}
