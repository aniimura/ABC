/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop48Cpt
import ABC3.Found.FrdI.Prop114
import ABC3.Found.FrdI.Prop44Gp
import ABC3.Found.FrdI.Prop44Phi

/-!
# [FrdI] Proposition 4.8, (iii), (b) —— `Frobenius-compact` の第 3 条

## ★★★★★★以前の見立ての誤り

台帳は「`σ` が恒等 ⟹ `Div^gp` の核に入る ⟹ **核が torsion であること**が要る」と
書いていたが、これは**遠回り**であった(`frobeniusCompact_cond3_of_eq` を見よ)。

★仮定は `∀ u ∈ 𝒪^×(A), ∃ k, (endConj θ u)^(d*k) = u^(c*k)` であり、
結論は `∀ u ∈ 𝒪^×(A), ∃ k, (endConj θ u)^k = u^k` である。
★★**`c = d` なら仮定の式そのものが結論**である。
★★★したがって `𝒪^×(A)` の torsion 性は**まったく関係が無い**。

## ★残る問い —— `c = d`

`Div^gp` を当てると `(d*k) • T(x) = (c*k) • x`(`T := Φ^gp(Base θ.inv)`、`x := Div^gp u`)。
★`T = id` が言えれば `(d*k) • x = (c*k) • x` となり、
`x` が torsion でなければ `c = d` が出る。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★1. `non-dilating` から「恒等」を出す —— 一般のモノイド版 -/

section MonoidGen

variable {M : Type w} [AddCommMonoid M]

/-- ★★★★**`isDivIdentity_of_forall_mprec` の一般形** ——
`Frobenioid` を経由しない、素の加法モノイド準同型についての主張。

★元の証明と同じ 2 段: `charMap` を `non-dilating` で恒等にし、
`toChar` の単射性(sharp)で `M` へ戻す。 -/
theorem addMonoidHom_eq_id_of_forall_mprec (hsharp : IsSharp M) (f : M →+ M)
    (hnd : IsNonDilating f) (h : ∀ x : M, x ≠ 0 → MPrec (f x) x) :
    f = AddMonoidHom.id M := by
  have hchar : charMap f = AddMonoidHom.id (MChar M) := by
    refine hnd (fun a ha => ?_)
    obtain ⟨x, rfl⟩ := toChar_surjective _ a
    have hx : x ≠ 0 := fun hz => ha.1 (by rw [hz, map_zero])
    rw [charMap_toChar]
    exact MPrec.map toChar (h x hx)
  refine AddMonoidHom.ext (fun x => ?_)
  refine toChar_injective_of_isSharp hsharp ?_
  have := congrArg (fun g : MChar M →+ MChar M => g (toChar x)) hchar
  rw [charMap_toChar] at this
  exact this

def addMonoidHom_eq_id_of_forall_mprec.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 41,
    item := "Proposition 1.14, (v) — non-dilating から恒等（一般形）",
    sectionId := "frdi-prop-1-14" }

end MonoidGen

/-! ## ★2. `Div^gp` と共役 -/

section BiratConj

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {G : Frobenioid P}

/-- ★`𝒞^birat` では同型の Frobenius 次数は 1。 -/
theorem biratDeg_of_isIso {A B : BiratCat P G} (f : A ⟶ B) [IsIso f] :
    ((biratDeg f : ℕ+) : ℕ) = 1 := by
  have h : (biratPre P G).degFr f = 1 := degFr_of_isIso (biratPre P G) f
  rw [biratPre_degFr] at h
  rw [h]
  rfl

/-- ★★★★★**共役の `Div^gp`** —— `Div^gp (endConj θ u) = Φ^gp(Base θ.inv) (Div^gp u)`。

★`θ.inv ≫ θ.hom = 𝟙` から `Div^gp θ.inv` が消える。 -/
theorem biratDivGp_endConj {A : BiratCat P G} (θ : A ≅ A) {u : End A}
    (hu : u ∈ OTimes (biratPre P G) A) :
    biratDivGp ((endConj θ u : End A) : A ⟶ A)
      = gpMap _ (Φ.map (biratBase θ.inv)) (biratDivGp ((u : A ⟶ A))) := by
  -- ★`u ≫ θ.hom`
  have h1 : biratDivGp ((u : A ⟶ A) ≫ θ.hom)
      = biratDivGp θ.hom + biratDivGp ((u : A ⟶ A)) := by
    have h := biratDivGp_comp' ((u : A ⟶ A)) θ.hom
    rw [gpMap_biratBase_of_baseIdentity hu.1.1, biratDeg_of_isIso θ.hom, one_smul] at h
    exact h
  -- ★`θ.inv ≫ (u ≫ θ.hom)`
  have h2 : biratDivGp (θ.inv ≫ ((u : A ⟶ A) ≫ θ.hom))
      = gpMap _ (Φ.map (biratBase θ.inv)) (biratDivGp ((u : A ⟶ A) ≫ θ.hom))
        + biratDivGp θ.inv := by
    have h := biratDivGp_comp' θ.inv ((u : A ⟶ A) ≫ θ.hom)
    have hdeg : ((biratDeg ((u : A ⟶ A) ≫ θ.hom) : ℕ+) : ℕ) = 1 := by
      have hiso : IsIso ((u : A ⟶ A) ≫ θ.hom) := by
        haveI : IsIso (u : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso u).mp hu.2
        infer_instance
      exact biratDeg_of_isIso _
    rw [hdeg, one_smul] at h
    exact h
  -- ★`𝟙 = θ.inv ≫ θ.hom` から `Div^gp θ.inv` を消す
  have h3 : gpMap _ (Φ.map (biratBase θ.inv)) (biratDivGp θ.hom) + biratDivGp θ.inv = 0 := by
    have h := biratDivGp_comp' θ.inv θ.hom
    rw [biratDeg_of_isIso θ.hom, one_smul, θ.inv_hom_id, biratDivGp_id] at h
    exact h.symm
  show biratDivGp (θ.inv ≫ ((u : A ⟶ A) ≫ θ.hom)) = _
  rw [h2, h1, map_add]
  rw [add_assoc, add_comm (gpMap _ (Φ.map (biratBase θ.inv)) (biratDivGp ((u : A ⟶ A)))),
    ← add_assoc, h3, zero_add]

def biratDivGp_endConj.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iii) — 共役の Div^gp",
    sectionId := "frdi-prop-4-4" }

/-- ★`𝒪^×` の冪の `Div^gp` は倍数。 -/
theorem biratDivGp_pow {A : BiratCat P G} {v : End A}
    (hv : v ∈ OTimes (biratPre P G) A) (n : ℕ) :
    biratDivGp (((v ^ n : End A)) : A ⟶ A) = n • biratDivGp ((v : A ⟶ A)) := by
  induction n with
  | zero =>
    show biratDivGp ((1 : End A) : A ⟶ A) = (0 : ℕ) • _
    rw [zero_nsmul]
    exact biratDivGp_id A
  | succ m ih =>
    have hvm : (v ^ m) ∈ OTimes (biratPre P G) A := pow_mem hv m
    show biratDivGp (((v ^ (m + 1) : End A)) : A ⟶ A) = _
    rw [pow_succ, biratDivGp_mul_otimes hvm hv, ih, succ_nsmul]

/-! ## ★3. `c = d` -/

/-- ★群では `m • y = n • y`(`m ≤ n`)なら `(n - m) • y = 0`。 -/
theorem nsmul_sub_eq_zero_of_eq {A : Type*} [AddCommGroup A] {y : A} {m n : ℕ}
    (hle : m ≤ n) (heq : m • y = n • y) : (n - m) • y = 0 := by
  have hsplit : (n - m) • y + m • y = n • y := by
    rw [← add_nsmul, Nat.sub_add_cancel hle]
  have h : (n - m) • y + m • y = 0 + m • y := by rw [hsplit, ← heq, zero_add]
  exact add_right_cancel h

/-- ★★★★★★**`c = d`** —— `Frobenius-compact` の第 3 条の本体。

★段は 3 つ:
1. 仮定に `Div^gp` を当て、`𝒪^×(A^birat) ↠ Φ^gp`(`Proposition 4.4, (iii)`)で
   **すべての** `x` について `(d*k) • σ(x) = (c*k) • x` を得る
2. `d*k ≥ 1` から `MPrec (σ x) x` が出るので、
   `non-dilating` で **`σ = id`**(`addMonoidHom_eq_id_of_forall_mprec`)
3. torsion でない `x₀` を 1 つ取れば `d*k = c*k`、したがって `c = d` -/
theorem birat_cd_eq
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (A : BiratCat P G) (θ : A ≅ A) (c d : ℕ+)
    (hnd : IsNonDilating (Φ.map (biratBase θ.inv)))
    (x₀ : Φ.val (P.toElem.obj (biratDown P G A)).base)
    (hx₀ : ∀ n : ℕ, 0 < n →
      n • toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x₀ ≠ 0)
    (hyp : ∀ u : End A, u ∈ OTimes (biratPre P G) A → ∃ k : ℕ+,
      ((endConj θ u) ^ (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)
        = (u ^ (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) : End A)) :
    c = d := by
  classical
  -- ★段 1
  have hkey : ∀ x : Φ.val (P.toElem.obj (biratDown P G A)).base,
      ∃ k : ℕ+, (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
          • gpMap _ (Φ.map (biratBase θ.inv))
              (toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x)
        = (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
          • toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x := by
    intro x
    have hmem : toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x ∈ phiBiratAt P G A := by
      rw [phiBiratAt_eq_top_of_divSurj P G hdivS A]
      exact AddSubgroup.mem_top _
    obtain ⟨u, hu, hux⟩ := hmem
    obtain ⟨k, hk⟩ := hyp u hu
    refine ⟨k, ?_⟩
    have hcj : (endConj θ u) ∈ OTimes (biratPre P G) A :=
      endConj_mem_otimes (biratPre P G) θ hu
    have h := congrArg (fun t : End A => biratDivGp ((t : A ⟶ A))) hk
    rw [biratDivGp_pow hcj, biratDivGp_pow hu, biratDivGp_endConj θ hu, hux] at h
    exact h
  -- ★段 2
  have hsharp : IsSharp (Φ.val (P.toElem.obj (biratDown P G A)).base) :=
    (P.divisorial _).2
  have hint : IsIntegralMonoid (Φ.val (P.toElem.obj (biratDown P G A)).base) :=
    (P.divisorial _).1.1
  have hσ : Φ.map (biratBase θ.inv) = AddMonoidHom.id _ := by
    refine addMonoidHom_eq_id_of_forall_mprec hsharp _ hnd (fun x _ => ?_)
    obtain ⟨k, hk⟩ := hkey x
    have hmono : (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) • (Φ.map (biratBase θ.inv) x)
        = (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) • x := by
      refine hint ?_
      rw [toGp_nsmul, toGp_nsmul, ← gpMap_toGp _ (Φ.map (biratBase θ.inv)) x]
      exact hk
    refine ⟨((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ), Nat.mul_pos c.2 k.2, ?_⟩
    refine ⟨(((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) - 1) • (Φ.map (biratBase θ.inv) x), ?_⟩
    have hpos : 0 < ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) := Nat.mul_pos d.2 k.2
    have h1 : Φ.map (biratBase θ.inv) x
          + (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) - 1) • (Φ.map (biratBase θ.inv) x)
        = (1 + (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) - 1)) • (Φ.map (biratBase θ.inv) x) := by
      rw [add_nsmul, one_nsmul]
    rw [h1, show 1 + (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) - 1)
        = ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) from by omega]
    exact hmono
  -- ★段 3
  obtain ⟨k, hk⟩ := hkey x₀
  rw [hσ] at hk
  have hk' : (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
        • toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x₀
      = (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
        • toGp (Φ.val (P.toElem.obj (biratDown P G A)).base) x₀ := by
    have hid : gpMap _ (AddMonoidHom.id (Φ.val (P.toElem.obj (biratDown P G A)).base))
        = AddMonoidHom.id _ := gpMap_id _
    rw [hid] at hk
    exact hk
  by_contra hne
  have hkpos : 0 < ((k : ℕ+) : ℕ) := k.2
  have hne' : ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) ≠ ((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) := by
    intro h
    exact hne (PNat.coe_injective (Nat.eq_of_mul_eq_mul_right hkpos h)).symm
  rcases Nat.lt_or_ge (((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ))
      (((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ)) with hlt | hge
  · exact hx₀ _ (by omega) (nsmul_sub_eq_zero_of_eq (le_of_lt hlt) hk')
  · have hgt : ((c : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) < ((d : ℕ+) : ℕ) * ((k : ℕ+) : ℕ) := by omega
    exact hx₀ _ (by omega) (nsmul_sub_eq_zero_of_eq (le_of_lt hgt) hk'.symm)

def birat_cd_eq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 23,
    item := "Definition 1.2, (iv) — Frobenius-compact の第 3 条: c = d",
    sectionId := "frdi-def-1-2-iv" }

end BiratConj

end ABC3.Found.FrdI
