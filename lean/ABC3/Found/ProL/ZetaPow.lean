/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.ProL.Decomposition

/-!
# ζ 乗写像 —— [FrdI] `Definition 2.8, (iii)` の本体

原文 (FrdI p.52):
> is a set-theoretic function, then we shall refer to as the map given by raising to the

原文 (FrdI p.53):
> power will always be bijective.]

## ★★測った要点(2026-08-19)

★原文の括弧書き「co-prime 型なら ζ 乗写像はつねに全単射」が本節の中身である。
成分ごとに見れば **pro-`l` 群で `l` と素な冪写像が全単射**であることに帰着する。

| 向き | 手筋 |
|---|---|
| 単射 | 各有限商は `l` 群 ⟹ 位数は `l` 冪、`n` と素 ⟹ 位数 1 |
| 全射 | 各有限商で `powCoprime`(在庫)、整合性は**単射性**から、`M` へ持ち上げ |

★★**位相的有限生成は要らない** —— 原文は `Definition 2.8, (ii)` でそれを仮定するが、
分解も ζ 乗の全単射性も**可換な副有限群なら成り立つ**。
★原文より弱い仮定で述べているので、原文の主張はこれの系である(逸脱の向きが安全側)。
-/

namespace ABC3.Found.ProL

open CategoryTheory

universe u

section ProLPow

variable {N : Type u} [CommGroup N] [TopologicalSpace N] [IsTopologicalGroup N]
  [CompactSpace N] [TotallyDisconnectedSpace N]

/-- ★pro-`l` 群では `l` と素な冪写像は**単射**。 -/
theorem pow_injective_of_proL {l n : ℕ} (hl : l.Prime)
    (hpro : ∀ U : OpenNormalSubgroup N, IsPGroup l (N ⧸ U.toSubgroup))
    (hcop : Nat.Coprime n l) : Function.Injective (fun x : N => x ^ n) := by
  intro x y hxy
  have h1 : (x * y⁻¹) ^ n = 1 := by
    have hx : x ^ n = y ^ n := hxy
    rw [mul_pow, inv_pow, hx, mul_inv_cancel]
  have h2 : x * y⁻¹ = 1 := by
    refine eq_of_forall_quotient_eq (asProfiniteGrp N) (fun U => ?_)
    obtain ⟨k, hk⟩ := hpro U (QuotientGroup.mk (x * y⁻¹))
    have hord1 : orderOf (QuotientGroup.mk (x * y⁻¹) : N ⧸ U.toSubgroup) ∣ n := by
      refine orderOf_dvd_of_pow_eq_one ?_
      rw [← QuotientGroup.mk_pow, h1]
      rfl
    have hord2 : orderOf (QuotientGroup.mk (x * y⁻¹) : N ⧸ U.toSubgroup) ∣ l ^ k :=
      orderOf_dvd_of_pow_eq_one hk
    have hcop2 : Nat.Coprime n (l ^ k) := Nat.Coprime.pow_right k hcop
    have hd1 : orderOf (QuotientGroup.mk (x * y⁻¹) : N ⧸ U.toSubgroup)
        ∣ Nat.gcd n (l ^ k) := Nat.dvd_gcd hord1 hord2
    rw [hcop2] at hd1
    have := orderOf_eq_one_iff.mp (Nat.dvd_one.mp hd1)
    rw [this]
    exact (QuotientGroup.mk_one _).symm
  exact mul_inv_eq_one.mp h2

/-- ★pro-`l` 群では `l` と素な冪写像は**全射**。 -/
theorem pow_surjective_of_proL {l n : ℕ} (hl : l.Prime)
    (hpro : ∀ U : OpenNormalSubgroup N, IsPGroup l (N ⧸ U.toSubgroup))
    (hcop : Nat.Coprime n l) : Function.Surjective (fun x : N => x ^ n) := by
  haveI : Fact l.Prime := ⟨hl⟩
  intro z
  -- ★各有限商では `powCoprime` が逆を与える
  have hcard : ∀ U : OpenNormalSubgroup N,
      Nat.Coprime (Nat.card (N ⧸ U.toSubgroup)) n := by
    intro U
    obtain ⟨k, hk⟩ := (IsPGroup.iff_card (p := l) (G := N ⧸ U.toSubgroup)).mp (hpro U)
    rw [hk]
    exact Nat.Coprime.pow_left k (Nat.Coprime.symm hcop)
  set f : ∀ U : OpenNormalSubgroup N, N ⧸ U.toSubgroup :=
    fun U => (powCoprime (n := n) (hcard U)).symm (QuotientGroup.mk z) with hf
  have hfn : ∀ U, (f U) ^ n = (QuotientGroup.mk z : N ⧸ U.toSubgroup) := by
    intro U
    exact (powCoprime (n := n) (hcard U)).apply_symm_apply (QuotientGroup.mk z)
  have hcompat : ∀ (U V : OpenNormalSubgroup N) (h : U.toSubgroup ≤ V.toSubgroup),
      QuotientGroup.map U.toSubgroup V.toSubgroup (MonoidHom.id N) h (f U) = f V := by
    intro U V h
    have hinj : Function.Injective (fun w : N ⧸ V.toSubgroup => w ^ n) := by
      intro a b hab
      obtain ⟨a', rfl⟩ := QuotientGroup.mk_surjective a
      obtain ⟨b', rfl⟩ := QuotientGroup.mk_surjective b
      exact (powCoprime (n := n) (hcard V)).injective hab
    refine hinj ?_
    show (QuotientGroup.map _ _ _ h (f U)) ^ n = (f V) ^ n
    rw [← map_pow, hfn U, hfn V]
    rfl
  obtain ⟨w, hw⟩ := exists_mk_eq_of_compatible (asProfiniteGrp N) f hcompat
  refine ⟨w, ?_⟩
  refine eq_of_forall_quotient_eq (asProfiniteGrp N) (fun U => ?_)
  show (QuotientGroup.mk (w ^ n) : N ⧸ U.toSubgroup) = QuotientGroup.mk z
  rw [QuotientGroup.mk_pow, hw U, hfn U]

/-- ★★★**pro-`l` 群では `l` と素な冪写像は全単射**。 -/
theorem pow_bijective_of_proL {l n : ℕ} (hl : l.Prime)
    (hpro : ∀ U : OpenNormalSubgroup N, IsPGroup l (N ⧸ U.toSubgroup))
    (hcop : Nat.Coprime n l) : Function.Bijective (fun x : N => x ^ n) :=
  ⟨pow_injective_of_proL hl hpro hcop, pow_surjective_of_proL hl hpro hcop⟩

end ProLPow

/-! ## ★2. ζ 乗写像 -/

section Zeta

variable {M : Type u} [CommGroup M] [TopologicalSpace M] [IsTopologicalGroup M]
  [CompactSpace M] [TotallyDisconnectedSpace M]

/-- ★★**[FrdI] Definition 2.8, (iii)** —— `ζ` が **co-prime 型**であること。 -/
def IsCoPrimeType (ζ : Nat.Primes → ℕ+) : Prop :=
  ∀ l : Nat.Primes, Nat.Coprime ((ζ l : ℕ)) l.1

variable (M) in
/-- ★★★**[FrdI] Definition 2.8, (iii)** —— **ζ 乗写像**。
`M[l]` の上で `ζ(l)` 乗する。 -/
noncomputable def zetaPow (ζ : Nat.Primes → ℕ+) (x : M) : M :=
  (decompEquiv M).symm (fun l => (decompEquiv M x l) ^ ((ζ l : ℕ)))

variable (M) in
/-- ★成分ごとの冪写像。 -/
noncomputable def zetaPowPi (ζ : Nat.Primes → ℕ+)
    (y : (l : Nat.Primes) → lPart M l.1) : (l : Nat.Primes) → lPart M l.1 :=
  fun l => (y l) ^ ((ζ l : ℕ))

variable (M) in
theorem zetaPow_eq (ζ : Nat.Primes → ℕ+) :
    zetaPow M ζ = (decompEquiv M).symm ∘ zetaPowPi M ζ ∘ (decompEquiv M) := rfl

variable (M) in
/-- ★★★★★**[FrdI] Definition 2.8, (iii) の括弧書き** ——
`ζ` が co-prime 型なら **ζ 乗写像は全単射**。 -/
theorem zetaPow_bijective (ζ : Nat.Primes → ℕ+) (hζ : IsCoPrimeType ζ) :
    Function.Bijective (zetaPow M ζ) := by
  have hcomp : ∀ l : Nat.Primes,
      Function.Bijective (fun x : lPart M l.1 => x ^ ((ζ l : ℕ))) := by
    intro l
    refine pow_bijective_of_proL (N := lPart M l.1) l.2 ?_ (hζ l)
    intro U
    exact isProL_lPartGrp (M := M) l.1 U
  have hpi : Function.Bijective (zetaPowPi M ζ) := by
    constructor
    · intro y z hyz
      funext l
      exact (hcomp l).1 (congrFun hyz l)
    · intro z
      refine ⟨fun l => ((hcomp l).2 (z l)).choose, ?_⟩
      funext l
      exact ((hcomp l).2 (z l)).choose_spec
  rw [zetaPow_eq]
  exact ((decompEquiv M).symm.bijective).comp (hpi.comp (decompEquiv M).bijective)

def zetaPow.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 52, item := "Definition 2.8, (iii) — ζ 乗写像",
    sectionId := "frdi-def-2-8" }

end Zeta

end ABC3.Found.ProL
