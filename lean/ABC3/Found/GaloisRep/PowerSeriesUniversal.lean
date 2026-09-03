/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.AdicEvalHom
import ABC3.Found.GaloisRep.TateDSeries
import ABC3.Meta.Claim
import Mathlib.RingTheory.AdicCompletion.Completeness

/-!
# 第 1127 ブロック —— **`PowerSeries` を万有な完備環として使う道具**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★これは何か

`Lemma 3.5` から `d + 1 < l` を外す道（§3 の枠の節点 3）の**降下の台**である。

☆第 1091 で測ったとおり、`R` を商体 `K` に移すと `q` の `I` 進収束が壊れる——
`K` は `I`-adic 完備でないからである。★しかし

    `PowerSeries A` は `(X)`-adic 完備（mathlib の instance）

であり、**係数環 `A` を何に取り替えても壊れない**。したがって

    `A₀ ≔ PowerSeries R`  →（係数を商体へ）→  `A₁ ≔ PowerSeries (FractionRing R)`

という道が通る。`A₁` では `l` も `1 − ζ^i` も可逆なので、
**既存の `hu` つきの補題がそのまま使える**。

☆最後に `A₀ → R`（`X ↦ q`、係数は恒等）で特殊化して `R` の中の等式に戻す。
★その特殊化が第 1126 の `evalAdicMapHom` である。

## ★★★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `map_mem_span_X_pow` | ★係数の写像は `(X)^n` を `(X)^n` に送る |
| `partialEvalMap_X` / `evalAdicMap_X` | ★`X` の値は `q` |
| `evalAdicMapHom_mem_pow` | ★★特殊化は `(X)^n` を `I^n` に送る（連続性） |
| `coeff_zero_of_mem_span_X` | ★`(X)` の元の定数項は `0` |
| `tateDXtail_mem` / `tateDYtail_mem` | ★`D` 側の尾も `I` に入る |
-/

namespace ABC3.Found.GaloisRep

open PowerSeries

/-! ## ★係数の写像は `(X)^n` を保つ -/

/-- ★★**係数環の写像は `(X)^n` を `(X)^n` に送る**——`PowerSeries.map` の連続性。 -/
theorem map_mem_span_X_pow {A B : Type} [CommRing A] [CommRing B] (ψ : A →+* B) (n : ℕ)
    {f : PowerSeries A} (hf : f ∈ (Ideal.span {(PowerSeries.X : PowerSeries A)}) ^ n) :
    PowerSeries.map ψ f ∈ (Ideal.span {(PowerSeries.X : PowerSeries B)}) ^ n := by
  rw [Ideal.span_singleton_pow, Ideal.mem_span_singleton] at hf ⊢
  obtain ⟨g, rfl⟩ := hf
  exact ⟨PowerSeries.map ψ g, by rw [map_mul, map_pow, PowerSeries.map_X]⟩

/-- ★`(X)` の元の定数項は `0`。 -/
theorem coeff_zero_of_mem_span_X {A : Type} [CommRing A] {f : PowerSeries A}
    (hf : f ∈ Ideal.span {(PowerSeries.X : PowerSeries A)}) :
    PowerSeries.coeff 0 f = 0 := by
  rw [Ideal.mem_span_singleton, PowerSeries.X_dvd_iff] at hf
  simpa using hf

/-! ## ★★★★特殊化 `A⟦q⟧ → R`（`X ↦ q`）の連続性 -/

section Specialize

variable {R S : Type} [CommRing R] [CommRing S] {I : Ideal R}

/-- ☆`X` の部分和は `n ≥ 2` でだけ `q` になる。 -/
theorem partialEvalMap_X (φ : S →+* R) (q : R) (n : ℕ) :
    partialEvalMap φ (PowerSeries.X : PowerSeries S) q n = if 1 < n then q else 0 := by
  simp [partialEvalMap, PowerSeries.coeff_X, Finset.sum_ite_eq']

/-- ★★**`X` の値は `q`**。 -/
theorem evalAdicMap_X [IsAdicComplete I R] (φ : S →+* R) (q : R) (hq : q ∈ I) :
    evalAdicMap φ (PowerSeries.X : PowerSeries S) q hq = q := by
  refine evalAdicMap_unique φ _ _ _ _ (fun n => ?_)
  rw [SModEq.sub_mem, partialEvalMap_X]
  by_cases h : 1 < n
  · simp [h]
  · rw [if_neg h, zero_sub]
    have hn : n ≤ 1 := by omega
    have hmem : q ∈ I ^ n := Ideal.pow_le_pow_right hn (by simpa using hq)
    simpa using neg_mem hmem

/-- ★★★★★**特殊化は連続である**——`(X)^n` を `I^n` に送る。

★これが `map_adicSum`・`map_tateXpairDF`（第 1125）が要求する仮説そのものである。 -/
theorem evalAdicMapHom_mem_pow [IsAdicComplete I R] (φ : S →+* R) (q : R) (hq : q ∈ I)
    (n : ℕ) (f : PowerSeries S)
    (hf : f ∈ (Ideal.span {(PowerSeries.X : PowerSeries S)}) ^ n) :
    evalAdicMapHom φ q hq f ∈ I ^ n := by
  rw [Ideal.span_singleton_pow, Ideal.mem_span_singleton] at hf
  obtain ⟨g, rfl⟩ := hf
  rw [map_mul, map_pow, evalAdicMapHom_apply, evalAdicMap_X]
  exact Ideal.mul_mem_right _ _ (Ideal.pow_mem_pow hq n)

end Specialize

/-! ## ★`D` 側の尾も `I` に入る -/

section Tail

variable {R : Type} [CommRing R] {I : Ideal R}

theorem tateDXtail_mem [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateDXtail u q hq ∈ I := by
  refine adicSum_mem _ _ ?_
  simpa using tateDXterm_mem_pow (I := I) (k := 1)
    (by simpa using Ideal.mul_mem_right u I hq)

theorem tateDYtail_mem [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateDYtail u q hq ∈ I := by
  refine adicSum_mem _ _ ?_
  simpa using tateDYterm_mem_pow (I := I) (k := 1)
    (by simpa using Ideal.mul_mem_right u I hq)

end Tail

/-! ## ★★★`PowerSeries` 側の単元 -/

section Units

variable {F : Type} [CommRing F]

/-- ★`C` は単元を単元に送る。 -/
theorem isUnit_C_of_isUnit {x : F} (hx : IsUnit x) : IsUnit (PowerSeries.C x : PowerSeries F) :=
  (PowerSeries.C (R := F)).isUnit_map hx

/-- ★★**体の上では `1 − C α` は `α ≠ 1` だけで単元**——`R` の側で要る `hu` が消える理由。 -/
theorem isUnit_one_sub_C {K : Type} [Field K] {α : K} (hα : α ≠ 1) :
    IsUnit (1 - PowerSeries.C α : PowerSeries K) := by
  have h : (1 : PowerSeries K) - PowerSeries.C α = PowerSeries.C (1 - α) := by
    rw [map_sub, map_one]
  rw [h]
  exact isUnit_C_of_isUnit (isUnit_iff_ne_zero.2 (sub_ne_zero.2 (Ne.symm hα)))

/-- ★★**体の上では `l` も単元**（`l ≠ 0` のとき）。 -/
theorem isUnit_natCast_powerSeries {K : Type} [Field K] {l : ℕ} (hl : (l : K) ≠ 0) :
    IsUnit ((l : ℕ) : PowerSeries K) := by
  have h : ((l : ℕ) : PowerSeries K) = PowerSeries.C ((l : ℕ) : K) := by simp
  rw [h]
  exact isUnit_C_of_isUnit (isUnit_iff_ne_zero.2 hl)

end Units

/-! ## ★出典の紐付け(`.src`) -/

def evalAdicMapHom_mem_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(q-展開の特殊化が連続であること——(X)^n を I^n に送る)",
    sectionId := "genell-def-3-3" }

def map_mem_span_X_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(係数環の取り替えは (X) 進位相を保つこと)",
    sectionId := "genell-def-3-3" }

def evalAdicMapHom_mem_pow.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsAdicComplete (span {X}) (PowerSeries R)(万有な完備環)"
      (.inMathlib "PowerSeries.instIsAdicCompleteSpanSingletonSetX") 1,
    .implicitStep
      ("★★**2026-09-01（第 1127）の設計**——第 1091 で「商体 `K` に移すと `q` の " ++
       "`I` 進収束が壊れる」と測った。☆`PowerSeries A` は `(X)`-adic 完備で、" ++
       "**係数環 `A` を何に取り替えても壊れない**——ここが抜け道である。" ++
       "★`A₀ = PowerSeries R → A₁ = PowerSeries (FractionRing R)` で `hu` を得て、" ++
       "`PowerSeries.map` の単射性で降ろし、`evalAdicMapHom` で `R` に特殊化する。") 4 ]

end ABC3.Found.GaloisRep
