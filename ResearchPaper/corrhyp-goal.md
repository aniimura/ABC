# CorrHyp 形式化ゴール

論文: `CorrHyp` = Shinichi Mochizuki, *Correspondences on Hyperbolic Curves*
([papers.json](papers.json) 登記、pdfPages 18、notationRisk **unmeasured** —— 逐語を使う前に必ず 260dpi で目視すること。まだ誰も目視していない)。
出典: `0_Source/Correspondences on Hyperbolic Curves.txt`(pdftotext 抽出、855 行、page marker 実測済み)。
「(comments)」版が併存するが本ゴールでは未参照。

## 節構成と対象範囲

本文は §0(Introduction)〜§6 の 7 節。§0 は `Theorem A`/`B`/`C` を持つが、
本文中に明記されている通りいずれも後続の定理の再掲である:

- `Theorem A`(p.1) ← 「cf. Lemma 4.1 and Theorem 4.2 in the text」
- `Theorem B`(p.2) ← 「Theorem B follows from Theorem 5.3 in the text」
- `Theorem C`(p.2) ← 「(given as Theorem 6.1 in the text)」

二重計上を避けるため §0 は対象に含めない(GenEll の §1 起点と同じ扱い)。
`Corollary` は本論文に 0 件。`Remark` は脚注的注記であり、GenEll の慣行
(Skeleton 化の対象は Definition/Proposition/Theorem/Lemma/Corollary)にならい対象外。

## ゴール(2026-09-04 設定、Skeleton 未着手なので進捗は全項目 0)

| § | 表題 | 項目(開始ページ) | 件数 | 進捗 |
|---|---|---|---|---|
| §1 | Basic Definitions | Def 1.1(p.3)相関 / Def 1.2(p.4)自明な相関 / Def 1.3(p.4)転置相関 / Def 1.4(p.4)相関の合成 / Def 1.5(p.4)双曲曲線の isogeny | 5 | 0/5 |
| §2 | Review of Results of Margulis and Takeuchi | Def 2.1(p.5)無限個の相関 / Def 2.2(p.5)Margulis-arithmetic / Def 2.3(p.5)Shimura-arithmetic / Prop 2.4(p.6)両arithmeticの同値 / Thm 2.5(p.7、[Marg]引用)arithmetic⇔無限個の相関 / Thm 2.6(p.7、[Take]引用)与えた型のarithmetic Xは有限個 | 6 | 0/6 |
| §3 | The Non-arithmetic Case | Def 3.1(p.8)hyperbolic core / Prop 3.2(p.8)Γ_Z⊆Comm(Γ_X) / Thm 3.3(p.8)非arithmeticなXにisogenousな曲線は有限個 | 3 | 0/3 |
| §4 | The Main Theorem | Lemma 4.1(p.9)isogenyの体拡大からの降下 / Thm 4.2(p.10)第一主定理(有限性、Theorem Aの本体) | 2 | 0/2 |
| §5 | Isogenies of General Curves | Lemma 5.1(p.11)hyperbolic coreのk上への降下 / Def 5.2(p.12)hyperbolic core Yの定義 / Thm 5.3(p.12)一般の曲線は自身のhyperbolic coreに一致(Theorem Bの本体) / Lemma 5.4(p.13)e_Yの下界 / Lemma 5.5(p.13)d,g_Y,r_Y,Σ_Yの有限性 / Lemma 5.6(p.13)非自明な軌跡の可構成性 / Thm 5.7(p.15)例外型でのhyperbolic coreの記述 | 7 | 0/7 |
| §6 | Interpretation of a Theorem of Royden | Thm 6.1(p.17、Roydenの定理からの帰結)M_{g,r}の非自明な自己同型・相関の非存在(Theorem Cの本体) | 1 | 0/1 |
| **合計** | | | **24** | **0/24** |

```
CorrHyp §1 0/5, §2 0/6, §3 0/3, §4 0/2, §5 0/7, §6 0/1 (合計 0/24)
```

## 次の一手

[[leaf-first-with-graph-feedback]] に従い、まず §1(Def 1.1-1.5、5 件、依存 0 —— 相関・isogeny の語彙そのもの)から
Skeleton を立てるのが妥当。§1 は §2 以降すべての前提であり層番号が最小。

着手前の注意:
- notationRisk が unmeasured なので、Skeleton 化に入る前に該当ページ(まず p.3-4)の 260dpi 目視を行うこと
  (`papers.json` の警告どおり。GenEll では行列の出現順が入れ替わる実害が出た——この論文固有の壊れ方はまだ不明)。
- G1(出典)・G4(axiom 禁止)を満たすこと。
- §2 の Theorem 2.5 / 2.6 は他論文([Marg]・[Take])からの引用結果で、本プロジェクトでは証明せず
  posit すべき対象になる可能性が高い(要検討——Interface 行き)。
