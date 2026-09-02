/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DiscIdentitySemistable
import ABC3.Meta.Claim

/-!
# 第 1393 ブロック —— **`v_p(Δ)` と `minDeltaExp` は 12 を法として等しい**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`p ∣ l` の場合に効く

`semistableAt_veluQuotientFull` に残る `p ∣ l`（良い素点）の場合では、
核の点が形式群に入りうるので `N` も `E′` も**整とは限らない**。

★しかし**どのモデルでも `v_p(Δ)` は `minDeltaExp` と 12 を法として等しい**
——モデルどうしは変数変換で結ばれ、`Δ ↦ u⁻¹²Δ` だからである。

☆したがって恒等式 `Δ(E)^l = Δ(E′)·N⁴` と `v_p(Δ(E)) = 0` から

    minDeltaExp p E′ ≡ v_p(Δ(E′)) = −4·v_p(N)   (mod 12)

であり、**`3 ∣ v_p(N)` さえ言えれば `12 ∣ minDeltaExp p E′`** が出る。

★`3 ∣ v_p(N)` の中身は「形式群の点では `v(2y + a₁x + a₃) = −3k`、
それ以外の奇位数の点では `= 0`」であり、そこが形式群の要る一点である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField ABC3.Meta

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★★★★★★★
**`v_p(Δ)` と `minDeltaExp` は 12 を法として等しい**——★**無条件**（第 1393）。

☆モデルどうしは変数変換で結ばれ、`Δ ↦ u⁻¹²Δ` だからである。 -/
theorem dvd_valAdd_Delta_sub_minDeltaExp (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0) :
    (12 : ℤ) ∣ valAdd p (Units.mk0 W.Δ hΔ) - minDeltaExp p W := by
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) W

  have h := minDeltaExp_eq p W hΔ C hC
  have hu : (Units.mk0 ((C • W).Δ) (variableChange_Delta_ne_zero W hΔ C))
      = C.u⁻¹ ^ 12 * Units.mk0 W.Δ hΔ := by
    refine Units.ext ?_
    show (C • W).Δ = _
    rw [WeierstrassCurve.variableChange_Δ]
    push_cast
    simp
  rw [h, hu, valAdd_mul, valAdd_pow, valAdd_inv]
  exact ⟨valAdd p C.u, by ring⟩

/-- ★★★★★★★★★★★★★★★★
**恒等式から `minDeltaExp` の 12 を法とする値が読める**——★**無条件**（第 1393）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`v_p(Δ(E)) = 0` の素点では `minDeltaExp p E′ ≡ −4·v_p(N) (mod 12)` である。
★★★これが `p ∣ l` の場合の道であり、残るのは **`3 ∣ v_p(N)`** だけになる。 -/
theorem dvd_minDeltaExp_of_disc_pow_eq (p : HeightOneSpectrum (𝓞 L))
    (E E' : WeierstrassCurve L) (hΔ : E.Δ ≠ 0) (hΔ' : E'.Δ ≠ 0)
    (hgood : valAdd p (Units.mk0 E.Δ hΔ) = 0)
    {N : L} (hN : N ≠ 0) {l : ℕ} (hid : E.Δ ^ l = E'.Δ * N ^ 4)
    (hN3 : (3 : ℤ) ∣ valAdd p (Units.mk0 N hN)) :
    (12 : ℤ) ∣ minDeltaExp p E' := by
  -- ★段 1: 恒等式に `valAdd` を当てる
  have hu : (Units.mk0 E.Δ hΔ) ^ l = (Units.mk0 E'.Δ hΔ') * (Units.mk0 N hN) ^ 4 := by
    refine Units.ext ?_
    simpa using hid
  have h1 := congrArg (valAdd p) hu
  rw [valAdd_pow, valAdd_mul, valAdd_pow, hgood, mul_zero] at h1
  -- ★段 2: `v_p(Δ(E′)) = −4·v_p(N)` は 12 で割れる
  obtain ⟨m, hm⟩ := hN3
  have h2 : (12 : ℤ) ∣ valAdd p (Units.mk0 E'.Δ hΔ') := ⟨-m, by omega⟩
  -- ★段 3: `minDeltaExp` は 12 を法として同じ
  obtain ⟨k, hk⟩ := dvd_valAdd_Delta_sub_minDeltaExp p E' hΔ'
  obtain ⟨n, hn⟩ := h2
  exact ⟨n - k, by omega⟩

/-! ## ★出典の紐付け(`.src`) -/

def dvd_valAdd_Delta_sub_minDeltaExp.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(Δ) と minDeltaExp は 12 を法として等しい。★無条件)",
    sectionId := "genell-lemma-3-5" }

def dvd_valAdd_Delta_sub_minDeltaExp.needs : List ProofObligation := []

def dvd_minDeltaExp_of_disc_pow_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(恒等式から minDeltaExp の 12 を法とする値が読める。★無条件)",
    sectionId := "genell-lemma-3-5" }

def dvd_minDeltaExp_of_disc_pow_eq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "minDeltaExp_eq(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.minDeltaExp_eq") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1393）**——`p ∣ l`（良い素点）の場合の道である。" ++
       "☆核の点が形式群に入ると `N` も `E′` も整でないが、" ++
       "`v_p(Δ)` と `minDeltaExp` は 12 を法として等しいので " ++
       "`3 ∣ v_p(N)` さえ言えれば `12 ∣ minDeltaExp p E′` が出る。" ++
       "★`3 ∣ v_p(N)` の中身は「形式群の点では `v(2y+a₁x+a₃) = −3k`、" ++
       "それ以外の奇位数の点では `= 0`」であり、そこが形式群の要る一点である。") 17 ]

end ABC3.Found.GaloisRep
