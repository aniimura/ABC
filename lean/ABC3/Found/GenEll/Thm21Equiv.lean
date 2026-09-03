/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Interface.GenEll.Thm21Setup
import ABC3.Found.GenEll.BDClass
import ABC3.Meta.Claim

/-!
# 第 1446 ブロック —— **★★★★★★★★★★★★★★★★★★★★★★★★
`Theorem 2.1` の同値 `(i) ⟺ (ii)` を両向きとも証明した**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.11–p.13。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

## ★★★★★★★★これは何か

`Skeleton/GenEll/Section2.lean` の `theorem_2_1` は **`(i) ⟹ (ii)` だけ**を取っていた。
そこには「実質は `(ii) ⟹ (i)` の側であり、それは取れていない」と明記してあった。

★本ブロックは **`(ii) ⟹ (i)` も証明する**。
☆機構は第 1443（`Skeleton/GenEll/Section2Converse.lean` の節点 8）と同じで、
`Interface/GenEll/Thm21Setup.lean` の `Thm21Data.cover` 欄——
**段 B の幾何**「`X(ℚ̄)^{≤d}` の中の列は部分列を取れば
ある compactly bounded subset に丸ごと入る」——を仮定すれば、
組み立て自身は完全に証明できる。

| 向き | 原文の言い方 | 本ブロック |
|---|---|---|
| `(i) ⟹ (ii)` | 「immediate from the definitions」 | ✅ BD-類の制限、1 行 |
| `(ii) ⟹ (i)` | 原文 p.11–p.13 の 3 ページ | ✅ **`cover` 欄から背理法で出る** |

## ★★★★数に入れないこと

`Thm21Data.cover` は `Interface` の**欄**であり、まだ構成されていない。
★したがってこれは「`Theorem 2.1` が絶対的に証明された」ことを意味しない
——`Theorem 3.8` / `Corollary 4.3` / `Corollary 4.4` が
`EllModuliData` の欄を `Interface` で受けているのとまったく同じ状況である。
☆欄を埋めるのに要るのは noncritical Belyi 写像（一般曲線版、Riemann–Roch）だけであり、
コンパクト性の側（`Found/GenEll/Thm21Extract.lean` の `exists_bad_seq_tendsto`）は
既に `sorry` 0 で手元にある。

## ★向きに注意

`BDge α β` は逐語どおり `α(x) − β(x) ≤ C`、すなわち `α ≲ β` である。
★abc の不等式はこちらの向き（`ht_ω ≤ (1+ϵ)(log-diff + log-cond) + C`）であり、
`BDle` ではない（`Found/GenEll/BDClass.lean` の `bdle_ne_bdge` が
2 つが別物であることを反例で示している）。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta ABC3.Interface.GenEll

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] Theorem 2.1**(Compactly Bounded Subsets and the ABC Conjecture)——
★**同値の両向き**（第 1446）。

原文 (GenEll p.11):
> Theorem 2.1. (Compactly Bounded Subsets and the ABC Conjecture) Let Σ be a finite set of prime numbers.

☆(i) `ht_ω ≲ (1+ϵ)·(log-diff_X + log-cond_D)` が `X(ℚ̄)^{≤d}` の上で成り立つ。
☆(ii) 同じ不等式が、**どんな compactly bounded subset `K_V` に制限しても**成り立つ。

★★`(i) ⟹ (ii)` は「定義から直ちに」——同じ定数 `C` がそのまま使える。
★★`(ii) ⟹ (i)` は背理法である。`(i)` が偽なら
`n < ht_ω(x_n) − (1+ϵ)(…)` なる列 `x_n` が取れる。
`cover` 欄でその部分列を呑み込む `K` を取ると、`K` の上では `≤ C` なのに
`n` は `C` を超える——矛盾。

★★★★**`cover` は `Interface` の欄であり未構成である**。
本定理は「その欄を仮定すれば同値が出る」ことを言っている。 -/
theorem theorem_2_1 {P : Type} (T : Thm21Data P)
    (htOmega logDiff logCond : P → ℝ) (eps : ℝ) :
    -- (i) `X(ℚ̄)^{≤d}` の上での abc 不等式
    BDge (fun x : ↥T.degLe => htOmega x.1)
         (fun x : ↥T.degLe => (1 + eps) * (logDiff x.1 + logCond x.1))
    ↔
    -- (ii) どんな compactly bounded subset に制限しても成り立つ
    (∀ KV : Set P, T.CB KV →
      BDge (fun x : ↥(KV ∩ T.degLe) => htOmega x.1)
           (fun x : ↥(KV ∩ T.degLe) => (1 + eps) * (logDiff x.1 + logCond x.1))) := by
  constructor
  · -- ★(i) ⟹ (ii): 「immediate from the definitions」——同じ `C` を使い回す。
    rintro ⟨C, hC⟩ KV _
    exact ⟨C, fun x => hC ⟨x.1, x.2.2⟩⟩
  · -- ★★(ii) ⟹ (i): 原文 p.11–p.13。背理法。
    intro hii
    by_contra h
    rw [BDge] at h
    push Not at h
    -- ☆各 `n` に「`n` より悪い点」を選んで列を作る。
    choose u hu using fun n : ℕ => h (n : ℝ)
    -- ★段 B の幾何: その部分列を呑み込む compactly bounded subset `K` が取れる。
    obtain ⟨K, φ, hmono, hCB, hmem⟩ := T.cover (fun n => (u n).1) (fun n => (u n).2)
    -- ★(ii) は `K` の上で定数 `C` を与える。
    obtain ⟨C, hC⟩ := hii K hCB
    obtain ⟨m, hm⟩ := exists_nat_gt C
    have hle := hC ⟨(u (φ m)).1, ⟨hmem m, (u (φ m)).2⟩⟩
    have hgt := hu (φ m)
    have hmm : (m : ℝ) ≤ ((φ m : ℕ) : ℝ) := by exact_mod_cast hmono.le_apply
    -- ★★`φ m` 番目の点は `K` に居るので `≤ C`、しかし作り方から `> φ m ≥ m > C`。
    simp only at hle hgt
    linarith

/-! ## ★出典の紐付け(`.src`) -/

def theorem_2_1.src : Source :=
  { paper := "GenEll", pdfPage := 11, item := "Theorem 2.1",
    sectionId := "genell-thm-2-1" }

def theorem_2_1.needs : List ProofObligation :=
  [ .citation "[ABC3]" "Thm21Data.cover(段 B の幾何。★Interface の欄で未構成)"
      (.inProject "ABC3" "ABC3.Interface.GenEll.Thm21Data") 12,
    .citation "[ABC3]" "exists_bad_seq_tendsto(コンパクト性で部分列を取る、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_bad_seq_tendsto") 1,
    .citation "[ABC3]" "thm_2_1_stepA / thm_2_1_stepB(数値の鎖、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.thm_2_1_stepA") 1,
    .otherPaper "[NCBelyi]"
      "Theorem 2.5(Belyi Maps Noncritical at Prescribed Points)の一般曲線版——★cover 欄を埋めるのに要る唯一のもの。ℙ¹ 版 Found/NCBelyi/Thm25P1.lean は sorry 0" 5,
    .implicitStep
      ("★★★★★**2026-09-02（第 1446）——`(ii) ⟹ (i)` も取れた**。" ++
       "☆以前は `Skeleton/GenEll/Section2.lean` の `theorem_2_1` が `(i) ⟹ (ii)` だけで、" ++
       "「実質は `(ii) ⟹ (i)` の側であり、それは取れていない」と書いてあった。" ++
       "★段 B の幾何を `Thm21Data.cover` 欄 1 本に絞ったところ、組み立ては完全に証明できた。" ++
       "★★ただし `cover` は未構成の欄である——" ++
       "`Theorem 3.8` / `Corollary 4.3` / `Corollary 4.4` が `EllModuliData` の欄を" ++
       "`Interface` で受けているのとまったく同じ状況であり、" ++
       "「絶対的に証明された」ことを意味しない。") 11 ]

end ABC3.Found.GenEll
