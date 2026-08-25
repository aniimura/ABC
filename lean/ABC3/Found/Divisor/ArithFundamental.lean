/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithFunctor

/-!
# 素点の基本等式（相対版）—— `Σ_{w | v} localDeg(w) = [M : L]`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

## ★★何のために要るか

`Theorem 6.4, (i)` の 5 圏のうち **`𝒞^rlf`** を算術で出すには
`Φ^rlf : MonoidOn 𝒟` の条件 (a)（節点 `rlf-flat`）が要る。

★★在庫の `isCharacteristicallyInjective_scMap_of_nsmul_retraction`
（`Def24ScTransport.lean`）により、これは

    h ∘ (Φ.map α) = k · id   なる h の構成

に帰着している。算術では **「上にある素点についての和」**

    h(y)(v) := Σ_{V | v} y(V)

を取れば `h(f(x))(v) = (Σ_{V | v} localDeg V) · x(v)` なので、
★**要るのは相対版の基本等式 `Σ_{V | v} localDeg(V) = [M : L]` ただ 1 本**である。

## ★★★中身（無限素点の側）

★mathlib の相対版（`InfinitePlace/Ramification.lean` の `card_isUnramified` ほか）は
**Galois 専用**で、一般の `L ⊆ M` 用の等式は無い（2026-08-25 に確認）。
そこで埋め込みの二重数え上げで作る:

    S := {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v}

| 数え方 | 結果 | 使う在庫 |
|---|---|---|
| `mk` のファイバー | `Σ_{w \| v} mult w` | `card_filter_mk_eq`（mathlib） |
| `ψ ↦ ψ ∘ algebraMap` のファイバー | `mult v · [M:L]` | `AlgHom.card`（mathlib） |

★★見積りは **150〜300 行**としていたが、mathlib に
`card_filter_mk_eq` / `comap_mk` / `AlgHom.card` が揃っていたので
**約 60 行**で済んだ。
-/

namespace ABC3.Found.Divisor

open NumberField InfinitePlace Finset

/-! ## ★1. 与えられた埋め込みの延長は `[M:L]` 個 -/

open scoped Classical in
/-- ★★`φ : L →+* ℂ` の `M` への延長はちょうど `[M : L]` 個。

★`ℂ` を `φ` で `L`-代数と見れば `M →ₐ[L] ℂ` そのものなので、
mathlib の `AlgHom.card` がそのまま効く。 -/
theorem card_filter_comp_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (φ : L →+* ℂ) :
    #{ψ : M →+* ℂ | ψ.comp (algebraMap L M) = φ} = Module.finrank L M := by
  classical
  letI : Algebra L ℂ := φ.toAlgebra
  have e : {ψ : M →+* ℂ // ψ.comp (algebraMap L M) = φ} ≃ (M →ₐ[L] ℂ) := by
    refine ⟨fun p => ⟨p.1, ?_⟩, fun f => ⟨f.toRingHom, ?_⟩, ?_, ?_⟩
    · intro x
      show p.1 (algebraMap L M x) = algebraMap L ℂ x
      rw [← RingHom.comp_apply, p.2]
      rfl
    · ext x
      show f (algebraMap L M x) = φ x
      exact f.commutes x
    · intro p; rfl
    · intro f; rfl
  rw [← Fintype.card_subtype, Fintype.card_congr e]
  exact AlgHom.card L M ℂ

/-! ## ★2. 二重数え上げ -/

open scoped Classical in
/-- ★★★★`v` の上にある埋め込みの個数は `mult v · [M:L]`。

★`ψ ↦ ψ ∘ algebraMap L M` のファイバーで数える。 -/
theorem card_filter_comap_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (v : InfinitePlace L) :
    #{ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v}
      = v.mult * Module.finrank L M := by
  classical
  rw [Finset.card_eq_sum_card_fiberwise
    (f := fun ψ : M →+* ℂ => ψ.comp (algebraMap L M))
    (t := {φ : L →+* ℂ | mk φ = v})
    (by
      intro ψ hψ
      rw [Finset.mem_coe, mem_filter_univ] at hψ
      rw [Finset.mem_coe, mem_filter_univ, ← comap_mk]
      exact hψ)]
  have hinner : ∀ φ ∈ ({φ : L →+* ℂ | mk φ = v} : Finset (L →+* ℂ)),
      #({a ∈ {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v} |
          a.comp (algebraMap L M) = φ}) = (fun _ : L →+* ℂ => Module.finrank L M) φ := by
    intro φ hφ
    rw [mem_filter_univ] at hφ
    have hset : ({a ∈ {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v} |
        a.comp (algebraMap L M) = φ} : Finset (M →+* ℂ))
        = ({ψ : M →+* ℂ | ψ.comp (algebraMap L M) = φ} : Finset (M →+* ℂ)) := by
      ext ψ
      simp only [mem_filter, mem_univ, true_and]
      constructor
      · rintro ⟨-, h2⟩; exact h2
      · intro h2
        refine ⟨?_, h2⟩
        rw [comap_mk, h2]; exact hφ
    show #({a ∈ {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v} |
        a.comp (algebraMap L M) = φ} : Finset (M →+* ℂ)) = Module.finrank L M
    rw [hset]
    exact card_filter_comp_eq L M φ
  rw [Finset.sum_congr rfl hinner, Finset.sum_const, smul_eq_mul, card_filter_mk_eq]

open scoped Classical in
/-- ★★★★★★**無限素点の相対基本等式** —— `Σ_{w | v} mult w = mult v · [M:L]`。

★同じ集合を `mk` のファイバーで数え直すだけである。 -/
theorem sum_mult_comap_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (v : InfinitePlace L) :
    ∑ w ∈ {w : InfinitePlace M | w.comap (algebraMap L M) = v}, w.mult
      = v.mult * Module.finrank L M := by
  classical
  rw [← card_filter_comap_eq L M v,
    Finset.card_eq_sum_card_fiberwise
      (f := fun ψ : M →+* ℂ => mk ψ)
      (t := {w : InfinitePlace M | w.comap (algebraMap L M) = v})
      (by
        intro ψ hψ
        rw [Finset.mem_coe, mem_filter_univ] at hψ
        rw [Finset.mem_coe, mem_filter_univ]
        exact hψ)]
  refine (Finset.sum_congr rfl ?_).symm
  intro w hw
  rw [mem_filter_univ] at hw
  have hset : ({a ∈ {ψ : M →+* ℂ | (mk ψ).comap (algebraMap L M) = v} | mk a = w}
      : Finset (M →+* ℂ)) = ({ψ : M →+* ℂ | mk ψ = w} : Finset (M →+* ℂ)) := by
    ext ψ
    simp only [mem_filter, mem_univ, true_and]
    constructor
    · rintro ⟨-, h2⟩; exact h2
    · intro h2
      exact ⟨by rw [h2]; exact hw, h2⟩
  rw [hset, card_filter_mk_eq]

/-! ## ★3. `localDeg` の形 -/

open scoped Classical in
/-- ★★★★★★**無限素点の相対基本等式**（`localDeg` の形）——
`Σ_{W | v} localDeg(W) = [M : L]`。

★`localDeg W = mult W / mult v` で `mult v ∣ mult W`（`mult_resInf_dvd`）なので、
上の等式を `mult v` で割ればよい。 -/
theorem sum_localDeg_inf_eq (L M : Type) [Field L] [Field M] [NumberField L] [NumberField M]
    [Algebra L M] (v : InfinitePlace L) :
    ∑ W ∈ {W : InfinitePlace M | resInf (L := L) W = v},
        localDeg (L := L) (Sum.inr W) = Module.finrank L M := by
  classical
  have hkey := sum_mult_comap_eq L M v
  have hmul : (∑ W ∈ {W : InfinitePlace M | resInf (L := L) W = v},
      localDeg (L := L) (Sum.inr W)) * v.mult
      = ∑ W ∈ {W : InfinitePlace M | resInf (L := L) W = v}, W.mult := by
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl (fun W hW => ?_)
    rw [mem_filter_univ] at hW
    have h1 : localDeg (L := L) (Sum.inr W) * (resInf (L := L) W).mult = W.mult :=
      localDeg_mul_mult W
    rw [← hW]
    exact h1
  have hset : ({W : InfinitePlace M | resInf (L := L) W = v} : Finset (InfinitePlace M))
      = ({W : InfinitePlace M | W.comap (algebraMap L M) = v} : Finset (InfinitePlace M)) := rfl
  rw [hset] at hmul
  rw [hkey] at hmul
  have hv : 0 < v.mult := by
    rcases mult_eq_one_or_two v with h | h <;> omega
  rw [hset]
  exact Nat.eq_of_mul_eq_mul_right hv (hmul.trans (by ring))

/-! ### ★出典の紐付け -/

def sum_mult_comap_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 無限素点の相対基本等式(Σ_{w|v} mult w = mult v · [M:L])",
    sectionId := "frdi-thm-6-4" }

def sum_localDeg_inf_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 無限素点の相対基本等式(localDeg の形)",
    sectionId := "frdi-thm-6-4" }

def sum_localDeg_inf_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "AlgHom.card(#(M →ₐ[L] ℂ) = [M:L])"
      (.inMathlib "AlgHom.card") 114,
    .citation "[mathlib]" "InfinitePlace.card_filter_mk_eq(#{φ | mk φ = w} = mult w)"
      (.inMathlib "NumberField.InfinitePlace.card_filter_mk_eq") 114,
    .citation "[mathlib]" "InfinitePlace.comap_mk((mk ψ).comap f = mk (ψ ∘ f))"
      (.inMathlib "NumberField.InfinitePlace.comap_mk") 114,
    .implicitStep
      ("★mathlib の相対版(card_isUnramified ほか)は Galois 専用なので、" ++
       "埋め込みの二重数え上げで作った") 114 ]

end ABC3.Found.Divisor
