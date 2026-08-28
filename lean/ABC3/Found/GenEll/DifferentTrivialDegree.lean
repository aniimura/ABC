/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LogDiffDescent
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★`base` を仮定から外す —— `[L:K] = 1` なら `𝔡 = ⊤`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

## ★★★★★★★★★★★★★★これは何か —— 測定のやり直しが当たった

`§9-901`（塔の帰納）は `base`（`[L:K] = 1` ⟹ `𝔡_{K,L} = ⊤`）を
**仮定として受けていた**。理由は当時の測定:

> mathlib に `differentIdeal … = ⊤` の判定が見当たらなかった
> （`Algebra.isUnramified_iff_differentIdeal_eq_top`・`differentIdeal_eq_top_iff`・
> `differentIdeal_self` いずれも無し、2026-08-28 実測）

★★★★**しかし `𝔡` を直接触らずに済む道があった**——**判別式**である。

## ★★★機構 —— 判別式で挟む

| 道具 | 役割 |
|---|---|
| `Algebra.finrank_eq_one_iff_bijective_algebraMap` | `[L:K] = 1` ⟺ `algebraMap K L` が全単射 |
| `NumberField.discr_eq_discr_of_algEquiv` | `K ≃ₐ[ℚ] L` ⟹ `disc K = disc L` |
| `natAbs_discr_eq_absNorm_differentIdeal_mul_natAbs_discr_pow` | `\|disc L\| = N(𝔡)·\|disc K\|^{[L:K]}` |
| `Ideal.absNorm_eq_one_iff` | `N(I) = 1` ⟺ `I = ⊤` |

★`[L:K] = 1` を 3 番目に入れると `|disc L| = N(𝔡)·|disc K|`、
2 番目から `|disc K| = |disc L| ≠ 0`、よって `N(𝔡) = 1`、
4 番目で `𝔡 = ⊤`。

## ★★★★測定の記録（やり直し）

★★★★★**「mathlib に無い」は「その名前で無い」でしかなかった**。
`differentIdeal` を主語にした判定は確かに無いが、
**判別式を主語にすれば全部あった**（2026-08-28 実測）。

★★在庫を引くとき、**求めている命題の主語を変えて引き直す**——
これは `§9-897`（消費側が `§9-VerticalBound` に既にあった）と同じ形の当たりである。

## ★これで何が外れたか

`§9-901`・`§9-904` の `base` 仮定が**消える**。
★残る仮定は `step`（次数 `p` の Galois 拡大で `p^p ∈ 𝔡`）だけである。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★★★★★★★次数 1 の拡大の different は単位イデアル -/

/-- ★★★★★★★★★★★★★★**`[L:K] = 1` なら `𝔡_{K,L} = ⊤`**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★`§9-901` が仮定として受けていた `base` そのものである。
★★`differentIdeal` を主語にした判定は mathlib に無いが、
**判別式を主語にすれば全部ある**——それで出した。 -/
theorem differentIdeal_eq_top_of_finrank_eq_one (K L : Type) [Field K] [NumberField K]
    [Field L] [NumberField L] [Algebra K L] (h : Module.finrank K L = 1) :
    differentIdeal (𝓞 K) (𝓞 L) = ⊤ := by
  have hb : Function.Bijective (algebraMap K L) :=
    Algebra.finrank_eq_one_iff_bijective_algebraMap.mp h
  have he : K ≃ₐ[ℚ] L :=
    AlgEquiv.ofRingEquiv (f := RingEquiv.ofBijective (algebraMap K L) hb)
      (by intro q; simp [RingEquiv.ofBijective])
  have hdisc : NumberField.discr K = NumberField.discr L :=
    NumberField.discr_eq_discr_of_algEquiv K he
  have hd := NumberField.natAbs_discr_eq_absNorm_differentIdeal_mul_natAbs_discr_pow K (𝓞 K) L (𝓞 L)
  rw [h, pow_one, hdisc] at hd
  have hne : (NumberField.discr L).natAbs ≠ 0 :=
    Int.natAbs_ne_zero.2 (NumberField.discr_ne_zero L)
  have h1 : Ideal.absNorm (differentIdeal (𝓞 K) (𝓞 L)) = 1 := by
    have hpos : 0 < (NumberField.discr L).natAbs := Nat.pos_of_ne_zero hne
    nlinarith [hd, hpos]
  exact Ideal.absNorm_eq_one_iff.mp h1

/-! ## ★★★★★★★★★★★★★★★`base` を外した形 -/

/-- ★★★★★★★★★★★★★**塔の帰納（`base` 無し）**。

★`§9-901` の `pow_mem_differentIdeal_of_pow_finrank` から `base` を外したもの。 -/
theorem pow_mem_differentIdeal_of_pow_finrank' (p : ℕ) [Fact p.Prime]
    (step : ∀ (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = p →
      ((p : 𝓞 L)) ^ p ∈ differentIdeal (𝓞 K) (𝓞 L)) :
    ∀ (k : ℕ) (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = p ^ k →
      ((p : 𝓞 L)) ^ (k * p) ∈ differentIdeal (𝓞 K) (𝓞 L) :=
  pow_mem_differentIdeal_of_pow_finrank p
    (fun K L _ _ _ _ _ _ h => differentIdeal_eq_top_of_finrank_eq_one K L h) step

/-- ★★★★★★★★★★★★★★**elementary claim の `p`-群の側（`base` 無し・一様な形）**。

    `[L:K] = p^k ≤ d`（Galois） ⟹ `p^N ∈ 𝔡_{L/K}`   （`N ≔ p·log_p d` 以上）

★右辺の `N` は `K` にも `L` にも依らない。 -/
theorem pow_mem_uniform_of_finrank_le' (p d : ℕ) [Fact p.Prime]
    (step : ∀ (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = p →
      ((p : 𝓞 L)) ^ p ∈ differentIdeal (𝓞 K) (𝓞 L))
    (hd : d ≠ 0) (N : ℕ) (hN : Nat.log p d * p ≤ N)
    (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
    [Algebra K L] [IsGalois K L] (k : ℕ) (hk : Module.finrank K L = p ^ k)
    (hle : Module.finrank K L ≤ d) :
    ((p : 𝓞 L)) ^ N ∈ differentIdeal (𝓞 K) (𝓞 L) :=
  pow_mem_uniform_of_finrank_le p d
    (fun K L _ _ _ _ _ _ h => differentIdeal_eq_top_of_finrank_eq_one K L h)
    step hd N hN K L k hk hle

/-- ★★★★★★★★★★★★★★★**elementary claim の消費側（`base` 無し）**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

    `[M:K] = p^k ≤ d`（Galois）、`L ⊆ M` ⟹ `log-diff(L) − log-diff(K) ≤ N·log p`

★★**残る仮定は `step`（次数 `p` の Galois 拡大で `p^p ∈ 𝔡`）だけ**である
——原文の「it suffices to consider the case where [L : K] = p」そのもの。 -/
theorem logDiffOfField_sub_le_uniform_of_le' (p d : ℕ) [Fact p.Prime]
    (step : ∀ (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = p →
      ((p : 𝓞 L)) ^ p ∈ differentIdeal (𝓞 K) (𝓞 L))
    (hd : d ≠ 0) (N : ℕ) (hN : Nat.log p d * p ≤ N)
    (K L M : Type) [Field K] [NumberField K] [Field L] [NumberField L]
    [Field M] [NumberField M] [Algebra K M] [IsGalois K M] [Algebra L M]
    (k : ℕ) (hk : Module.finrank K M = p ^ k) (hle : Module.finrank K M ≤ d) :
    logDiffOfField L - logDiffOfField K ≤ (N : ℝ) * Real.log p :=
  logDiffOfField_sub_le_uniform_of_le p d
    (fun K L _ _ _ _ _ _ h => differentIdeal_eq_top_of_finrank_eq_one K L h)
    step hd N hN K L M k hk hle

/-! ## ★出典の紐付け(`.src`) -/

def differentIdeal_eq_top_of_finrank_eq_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——[L:K] = 1 なら 𝔡 = ⊤)",
    sectionId := "genell-prop-1-7" }

def pow_mem_differentIdeal_of_pow_finrank'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——塔の帰納、base 無し)",
    sectionId := "genell-prop-1-7" }

def pow_mem_uniform_of_finrank_le'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim の p-群の側——base 無し・一様な形)",
    sectionId := "genell-prop-1-7" }

def logDiffOfField_sub_le_uniform_of_le'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim の消費側——base 無し)",
    sectionId := "genell-prop-1-7" }

def differentIdeal_eq_top_of_finrank_eq_one.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Algebra.finrank_eq_one_iff_bijective_algebraMap"
      (.inMathlib "Algebra.finrank_eq_one_iff_bijective_algebraMap") 1,
    .citation "[mathlib]" "NumberField.discr_eq_discr_of_algEquiv"
      (.inMathlib "NumberField.discr_eq_discr_of_algEquiv") 1,
    .citation "[mathlib]"
      "natAbs_discr_eq_absNorm_differentIdeal_mul_natAbs_discr_pow(判別式と different)"
      (.inMathlib "NumberField.natAbs_discr_eq_absNorm_differentIdeal_mul_natAbs_discr_pow") 1,
    .citation "[mathlib]" "Ideal.absNorm_eq_one_iff(N(I) = 1 ⟺ I = ⊤)"
      (.inMathlib "Ideal.absNorm_eq_one_iff") 1,
    .implicitStep
      ("★★★★★測定のやり直し(2026-08-28): §9-901 は「mathlib に " ++
       "differentIdeal … = ⊤ の判定が無い」として base を仮定に置いたが、" ++
       "**それは「その名前で無い」でしかなかった**。" ++
       "differentIdeal を主語にした判定は確かに無いが、" ++
       "**判別式を主語にすれば全部あった**") 4,
    .implicitStep
      ("★★在庫を引くとき、求めている命題の**主語を変えて引き直す**" ++
       "——§9-897(消費側が §9-VerticalBound に既にあった)と同じ形の当たりである") 3,
    .implicitStep
      ("★これで §9-901・§9-904 の base 仮定が消えた。" ++
       "残る仮定は step(次数 p の Galois 拡大で p^p ∈ 𝔡)だけである" ++
       "——原文の「it suffices to consider the case where [L : K] = p」そのもの") 3 ]

end ABC3.Found.GenEll
