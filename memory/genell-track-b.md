---
name: genell-track-b
description: [GenEll] トラック(B)の長期ゴールと、実測で分かった律速要因・臨界パスの到達点。
metadata:
  type: project
---

**北極星(B トラック)**: `[GenEll] Corollary 4.4` まで、IUTchIV が要る範囲の [GenEll] を
`Found/` に `sorry` 無しで載せる。**同じ土台に立つ第二の頂は `Theorem 2.1`**
(「abc ⟺ ℙ¹∖{0,1,∞} 版 abc」。★**IUT を一切使わない**)。
計画本体は `D:\Math_ABC3\ResearchPaper\genell-goal.md`、進捗計数は
`node tools/genell-progress.mjs`(分母 = IUT 側からの需要の推移閉包 **24 件**)。

**★律速要因は 2 つの不在である(2026-08-16 実測、探索範囲つき)**:

| 塊 | 件数 | 律速 | mathlib |
|---|---:|---|---|
| §1(高さ) | 9 | Arakelov 理論(算術直線束) | **0 件** |
| §2(Thm 2.1) | 1 | §1 ＋ noncritical Belyi | Belyi **0 件**(★[Mzk1] は `0_Source` にある) |
| §3(Galois 作用) | 9 | l 捩れ点への Galois 表現 | **0 件**(2-torsion polynomial のみ) |
| §4(Cor 4.4) | 5 | §3 | — |

★**したがって B は「IUT の転写」ではなく、古典的数学を 2 本作る仕事**である。
`Mathlib/NumberTheory/Height/` は 2006 行あるが **ℙⁿ の Weil 高さ**であって
算術直線束の枠組みではない——「高さがあるから軽い」という着手前の見込みは**外れた**。

**Why:** 分母 24 のうち完成は 0 件。カウンタが動かないのは怠慢ではなく、
`.src` の 2 値規則(条つきは数えない)と上表の律速による。

**How to apply:** 着手するなら唯一カウンタを動かせるのは `Lemma 3.1`(残りは (iv))。
(iv) を 9 段に割り、**6 段が済んでいる**(`Found/GenEll/` の 4 ファイル、すべて sorry 0):

- 済: 柱1(`SL(2,ℤ_l)` が基本行列で生成。★局所環で成り立つので位相不要)/
  1a(合同核が加法群)/ 1b-有限(`ℤ/l^{n+1}` で `(l^n)² = 0`)/
  1b-共役(**共役作用が随伴作用に落ちる**)/ 3(`𝔰𝔩₂(𝔽_l)` の既約性。要るのは `2 ≠ 0` だけ)/
  4a(`Σ_{i<l}(1+D)^i = 0`。★**`l ≥ 5` はここの `C(l,3)` 1 箇所だけで効く**、負の対照つき)
- 未: (A) の帰納法の組み立て / 4b(`Ad(u)=1+D`・`D³=0` の同定)/ **(B) 逆極限(唯一の位相依存)**

関連: [[abc3-plan-two-track]] / [[challenger-audit-without-context]] / [[genell-bd-class-direction]]
