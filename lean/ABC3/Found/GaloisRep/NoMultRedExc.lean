/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.NorthcottHtJ
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★乗法還元を持たない類の例外集合は有限（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

## ★★★★★★★★★★★★★★★★★★これは何か —— `Lemma 3.7` の**第 2 の主張の核**

原文 `Lemma 3.7` の証明:

> if condition (b) is satisfied, and, moreover `E_L` has **no primes of multiplicative
> reduction** [so `deg∞([E_L]) = 0`], then the fact that `[E_L] ∈ K_V` implies that
> `ht∞([E_L])` is **bounded** [independently of `L`, `E_L`, `l`], so, by Proposition 1.4, (iv),
> `[E_L]` belongs to some [fixed] **finite** exceptional set `Exc_d`

★★**本ファイルはこれを証明する**——★**無条件**である。

`Interface/GenEll/EllModuli.lean` の `noMultRedExc`／`galoisFinite_noMultRedExc` 欄の
中身にあたる。

## ★機構 —— 2 つとも手元にある

| 段 | 内容 | 出どころ |
|---|---|---|
| 1 | `deg∞ = 0` なら `h(j)` の**有限素点側が `0`** | ★`htFinJ_le_degInfOf`（`§9-1003`、**無条件**）＋ `htFinJ_nonneg` |
| 2 | compactly bounded なら**アルキメデス側が有界** | ★`‖σ(j)‖ ≤ M` から `h∞(j) ≤ max(log M, 0)` |
| 3 | したがって `h(j)` が有界 ⟹ **有限**（Northcott） | ★`northcott_htJ_image`（`§9-1005`、**無条件**） |

★★★`ht∞ ≔ h(j)` の読み替え（`§9-1009` の逸脱）の下で、原文の
「`ht∞([E_L])` is bounded」がそのまま段 1＋2 になる。

## ★「Galois-finite」であること（`Example 1.3, (i)`）

原文 p.5:「If `E ⊆ X(ℚ̄)` is a subset such that each `E^{≤d}` [where `d` ranges over the
positive integers] is finite, then we shall say that `E` is **Galois-finite**」。

★本ファイルの主張は**次数 `d` を止めた族について有限**という形なので、
`d` を動かせばそのまま Galois-finite である。

## ☆`Lemma 3.7` に残るもの

☆第 3 の主張（`(a) ∨ (b)` かつ l-cyclic なら `[E] ∈ Exc`）には `lcyclicExc` が要り、
それは `Lemma 3.5` を経由する（残る唯一の入力は Faltings 高さの同種写像評価——
`Found/GaloisRep/Lemma35Concrete.lean` に型で固定してある）。
★第 1 の主張は `Found/GaloisRep/Lemma37A.lean` で**取れている**。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★段 1 —— 乗法還元が無ければ有限素点側は `0` -/

/-- ★★★★★★**`deg∞ = 0` なら `h(j)` の有限素点側は `0`**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★`htFinJ ≤ deg∞`（`§9-1003`、**無条件**）と `htFinJ ≥ 0` から挟む。 -/
theorem htFinJ_eq_zero_of_degInfOf_eq_zero (E : WeierstrassCurve L) [E.IsElliptic]
    (h : degInfOf L E = 0) : htFinJ L E = 0 :=
  le_antisymm (h ▸ htFinJ_le_degInfOf E) (htFinJ_nonneg E)

/-! ## ★★★★★★段 2 —— compactly bounded ならアルキメデス側が有界 -/

/-- ★★★★★★**`‖σ(j)‖ ≤ M` なら `h∞(j) ≤ max(log M, 0)`**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★これが原文の「the fact that `[E_L] ∈ K_V` implies that `ht∞([E_L])` is bounded」の
アルキメデス側である（compactly bounded subset は各アルキメデス素点で
コンパクト集合に入ることを課す）。
★★`j = 0` でも `Real.log 0 = 0` なので破綻しない。 -/
theorem htArchJ_le_of_bound (E : WeierstrassCurve L) [E.IsElliptic] (M : ℝ)
    (hM : ∀ σ : L →+* ℂ, ‖σ E.j‖ ≤ M) :
    htArchJ L E ≤ max (Real.log M) 0 := by
  have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
  have hterm : ∀ σ : L →+* ℂ, max (Real.log ‖σ E.j‖) 0 ≤ max (Real.log M) 0 := by
    intro σ
    rcases (norm_nonneg (σ E.j)).lt_or_eq with hx | hx
    · exact max_le_max (Real.log_le_log hx (hM σ)) le_rfl
    · rw [← hx]
      simp [Real.log_zero]
  have hsum : ∑ σ : (L →+* ℂ), max (Real.log ‖σ E.j‖) 0
      ≤ ∑ _σ : (L →+* ℂ), max (Real.log M) 0 := Finset.sum_le_sum (fun σ _ => hterm σ)
  rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, NumberField.Embeddings.card L ℂ] at hsum
  rw [htArchJ, div_le_iff₀ hd]
  linarith

/-! ## ★★★★★★★★★★★★★★★★★★段 3 —— 例外集合は有限 -/

/-- ★★★★★★★★★★★★★★★★★★**乗法還元を持たない compactly bounded な類は有限**。

原文 (GenEll p.18):
> Lemma 3.7. (Finite Exceptional Sets) Let

★★★**無条件**である。`Interface/GenEll/EllModuli.lean` の
`noMultRedExc`／`galoisFinite_noMultRedExc` 欄の中身にあたる。

★仮定は原文自身のもの:

| 仮定 | 原文 |
|---|---|
| `hnomult` | 「`E_L` has **no primes of multiplicative reduction** [so `deg∞([E_L]) = 0`]」 |
| `harch` | 「`[E_L] ∈ K_V`」（compactly bounded subset のアルキメデス側） |
| `hdeg` | `M_ell(ℚ̄)^{≤d}` |

★★次数 `d` を止めた族について有限なので、`d` を動かせば
`Example 1.3, (i)` の意味で **Galois-finite** である。 -/
theorem finite_j_of_noMultRed (d : ℕ) {P : Type}
    (fld : P → IntermediateField ℚ ℂ) (hnf : ∀ p, NumberField (fld p))
    (hdeg : ∀ p, Module.finrank ℚ (fld p) ≤ d)
    (E : ∀ p, WeierstrassCurve (fld p)) (hE : ∀ p, (E p).IsElliptic)
    (M : ℝ)
    (hnomult : ∀ p, haveI := hnf p; haveI := hE p; degInfOf (fld p) (E p) = 0)
    (harch : ∀ p, haveI := hnf p; haveI := hE p; ∀ σ : (fld p) →+* ℂ, ‖σ (E p).j‖ ≤ M) :
    ((fun p : P => haveI := hE p; (((E p).j : fld p) : ℂ)) '' Set.univ).Finite := by
  refine Set.Finite.subset
    (northcott_htJ_image d fld hnf hdeg E hE
      (fun p => haveI := hnf p; haveI := hE p; htJ (fld p) (E p)) (fun p => rfl)
      (max (Real.log M) 0)) ?_
  refine Set.image_mono (fun p _ => ?_)
  haveI := hnf p; haveI := hE p
  show htJ (fld p) (E p) ≤ max (Real.log M) 0
  rw [htJ, htFinJ_eq_zero_of_degInfOf_eq_zero (E p) (hnomult p), zero_add]
  exact htArchJ_le_of_bound (E p) M (harch p)

/-! ## ★出典の紐付け(`.src`) -/

def htFinJ_eq_zero_of_degInfOf_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(deg∞ = 0 なら h(j) の有限素点側は 0。★無条件)",
    sectionId := "genell-lemma-3-7" }

def finite_j_of_noMultRed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Lemma 3.7(第 2 の主張の核——乗法還元を持たない compactly bounded な類は有限。★無条件)",
    sectionId := "genell-lemma-3-7" }

def finite_j_of_noMultRed.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htFinJ_le_degInfOf(h_fin(j) ≤ deg∞——無条件、§9-1003)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.htFinJ_le_degInfOf") 3,
    .citation "[ABC3]" "northcott_htJ_image(htJ が有界な曲線の j は有限個、§9-1005)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.northcott_htJ_image") 4,
    .otherPaper "[GenEll]"
      ("Proposition 1.4, (iv)(Northcott)——★済(Found/GenEll/Prop14.lean、§1 で計上済み)。" ++
       "本ファイルは §9-1005 の northcott_htJ_image 経由で使っている") 3,
    .otherPaper "[GenEll]"
      ("Example 1.3, (ii)(compactly bounded subset)——★本ファイルはそのアルキメデス側の" ++
       "帰結だけを `harch : ∀ σ, ‖σ(j)‖ ≤ M` として受けている。" ++
       "☆compactly bounded の定義そのもの(非アルキメデス側の Kv も含む)は取っていない") 5,
    .implicitStep
      ("★★★★★★★★★★到達点(2026-08-29、第 567): Lemma 3.7 の**第 2 の主張の核が" ++
       "無条件で取れた**。★原文の「ht∞([E_L]) is bounded」は、ht∞ ≔ h(j) の読み替え" ++
       "(§9-1009 の逸脱)の下で『有限素点側が 0』＋『アルキメデス側が有界』になり、" ++
       "どちらも手元にある。★★これで Lemma 3.7 の 3 主張のうち**第 1・第 2 が取れた**" ++
       "(第 1 は Found/GaloisRep/Lemma37A.lean、§9-1013)") 9,
    .folklore
      ("☆第 3 の主張には lcyclicExc が要り、それは Lemma 3.5 を経由する。" ++
       "★Lemma 3.5 に残る唯一の入力は Faltings 高さの同種写像評価であり、" ++
       "Found/GaloisRep/Lemma35Concrete.lean に型で固定してある") 9 ]

end ABC3.Found.GaloisRep
