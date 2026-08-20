/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm51Pic

/-!
# [FrdI] Theorem 5.1, (ii) —— 射の存在条件

原文 (FrdI p.97):
> in C of Frobenius degree d such that Base(

原文 (FrdI p.97):
> [where, by abuse of notation, we denote by z|AD the image of

★証明は 2 段に分かれる:

1. **底変換**(`picPullBack`)—— `θ : A_𝒟 → A'_𝒟` に沿って `Pic` の対象を引き戻すと、
   類は `Φ(θ)` でちょうど写る。★pull-back 射の普遍性だけで出る。
2. **同底の場合**(`picHom_sameBase_of_mem` / `_mem_of`)——
   `Definition 1.3, (iv)` の分解 `φ = γ ≫ β₀ ≫ α`(Frobenius・pre-step・pull-back)に対し、
   `γ` が `d` 倍(`hasCls_frobType`)、`β₀` が `Div` ぶんのずれ(`hasCls_preStep`)、
   `α` は**底が同型な pull-back なので同型**(`isIso_of_pullBack_base_isIso`)。

★1 の前提として「**pull-back 射に沿って Frobenius-trivial 性が引き戻される**」
(`isFrobeniusTrivial_pullBack`)が要る。これは原文が書いていない一歩だが、
pull-back の普遍性で `ζ` を持ち上げるだけである。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section PhiTools

variable {D : Type u} [Category.{v} D] {Φ : MonoidOn.{v, u, w} D}

theorem phi_map_comp_inv {X Y : D} (g : X ⟶ Y) [IsIso g] (z : Φ.val X) :
    Φ.map g (Φ.map (inv g) z) = z := by
  have h : Φ.map g (Φ.map (inv g) z) = Φ.map (g ≫ inv g) z := (Φ.map_comp _ _ z).symm
  rw [h, IsIso.hom_inv_id]
  exact Φ.map_id _ z

end PhiTools

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. pull-back 射の普遍性 -/

/-- ★pull-back 射の普遍性(存在)。 -/
theorem pullBack_lift {A B : C} {π : A ⟶ B} (h : IsPullBack P π) {X : C} (f : X ⟶ B)
    (g : (P.toElem.obj X).base ⟶ (P.toElem.obj A).base)
    (hcomm : P.Base f = g ≫ P.Base π) :
    ∃ f' : X ⟶ A, f' ≫ π = f ∧ P.Base f' = g := by
  obtain ⟨f', hf'⟩ := (h X).2 ⟨(f, g), hcomm⟩
  have h1 := congrArg (fun t => t.1) hf'
  exact ⟨f', congrArg Prod.fst h1, congrArg Prod.snd h1⟩

/-- ★pull-back 射の普遍性(一意性)。 -/
theorem pullBack_uniq {A B : C} {π : A ⟶ B} (h : IsPullBack P π) {X : C} {f₁ f₂ : X ⟶ A}
    (h1 : f₁ ≫ π = f₂ ≫ π) (h2 : P.Base f₁ = P.Base f₂) : f₁ = f₂ :=
  (h X).1 (Subtype.ext (Prod.ext h1 h2))

/-- ★底が同型な pull-back 射は同型。 -/
theorem isIso_of_pullBack_base_isIso {A B : C} {α : A ⟶ B} (h : IsPullBack P α)
    (hb : IsIso (P.Base α)) : IsIso α := by
  haveI := hb
  obtain ⟨α', hα'1, hα'2⟩ := pullBack_lift h (𝟙 B) (inv (P.Base α)) (by
    rw [P.Base_id, IsIso.inv_hom_id])
  refine ⟨α', ?_, hα'1⟩
  refine pullBack_uniq h (X := A) (f₁ := α ≫ α') (f₂ := 𝟙 A) ?_ ?_
  · rw [Category.assoc, hα'1, Category.comp_id, Category.id_comp]
  · rw [P.Base_comp, hα'2, IsIso.hom_inv_id, P.Base_id]

/-- ★★`Definition 1.3, (i), (c)` の本質的全射性 —— 底の任意の射は pull-back 射で実現される。 -/
theorem exists_pullBack_over (G : Frobenioid P) (A' : C) {Y : D}
    (θ : Y ⟶ (P.toElem.obj A').base) :
    ∃ (A₀ : C) (π : A₀ ⟶ A') (ρ : (P.toElem.obj A₀).base ≅ Y),
      IsPullBack P π ∧ ρ.hom ≫ θ = P.Base π := by
  haveI := G.core.plBkEquiv A'
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage
    (F := plBkOverFunctor P A') (Over.mk θ)
  exact ⟨Z.left.obj, Z.hom.hom, (Over.forget _).mapIso e, Z.hom.property, Over.w e.hom⟩

/-! ## ★2. Frobenius-trivial 性は pull-back で引き戻される -/

/-- ★★★**pull-back 射に沿って Frobenius-trivial 性は引き戻される**。

★原文は `Theorem 5.1, (iv)` で「(iii) から直ちに」と書くが、
その「直ちに」の中身がこれである。★`ζ` を pull-back の普遍性で持ち上げるだけ。 -/
theorem isFrobeniusTrivial_pullBack (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A₀ A' : C} (π : A₀ ⟶ A') (hπ : IsPullBack P π) (hA' : IsFrobeniusTrivial P A') :
    IsFrobeniusTrivial P A₀ := by
  obtain ⟨ζ, hdeg, hprop⟩ := hA'
  have hπd : P.Div π = 0 := (G.core.pullBackLB π hπ).1.2
  have hπl : P.degFr π = 1 := (G.core.pullBackLB π hπ).2
  have key : ∀ n : ℕ+, ∃ κ : A₀ ⟶ A₀,
      κ ≫ π = π ≫ ((ζ n : End A') : A' ⟶ A') ∧ P.Base κ = P.Base (𝟙 A₀) := by
    intro n
    refine pullBack_lift hπ _ (P.Base (𝟙 A₀)) ?_
    rw [P.Base_comp, show P.Base ((ζ n : End A') : A' ⟶ A') = P.Base (𝟙 A') from (hprop n).1,
      P.Base_id, Category.comp_id, P.Base_id, Category.id_comp]
  choose κ hκc hκb using key
  have huniq : ∀ (n : ℕ+) (f : A₀ ⟶ A₀), f ≫ π = π ≫ ((ζ n : End A') : A' ⟶ A') →
      P.Base f = P.Base (𝟙 A₀) → f = κ n := by
    intro n f h1 h2
    exact pullBack_uniq hπ (h1.trans (hκc n).symm) (h2.trans (hκb n).symm)
  have hdegκ : ∀ n : ℕ+, P.degFr (κ n) = n := by
    intro n
    have h := congrArg P.degFr (hκc n)
    rw [P.degFr_comp, P.degFr_comp, hπl, one_mul, mul_one, hdeg n] at h
    exact h
  have hdivκ : ∀ n : ℕ+, P.Div (κ n) = 0 := by
    intro n
    have h := congrArg P.Div (hκc n)
    rw [P.Div_comp, P.Div_comp, hπd, map_zero, smul_zero, add_zero,
      show P.Div ((ζ n : End A') : A' ⟶ A') = 0 from (hprop n).2.1.2, map_zero, zero_add,
      hπl] at h
    show P.Div (κ n) = 0
    rw [← h]
    show _ = ((1 : ℕ+) : ℕ) • P.Div (κ n)
    rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  refine ⟨{ toFun := fun n => (κ n : End A₀), map_one' := ?_, map_mul' := ?_ }, hdegκ, ?_⟩
  · exact (huniq 1 (𝟙 A₀) (by rw [Category.id_comp, map_one]; exact (Category.comp_id _).symm)
      rfl).symm
  · intro m n
    have hmul : ((ζ (m * n) : End A') : A' ⟶ A')
        = ((ζ n : End A') : A' ⟶ A') ≫ ((ζ m : End A') : A' ⟶ A') := by
      rw [map_mul]; rfl
    show κ (m * n) = (κ n : A₀ ⟶ A₀) ≫ κ m
    refine (huniq (m * n) (κ n ≫ κ m) ?_ ?_).symm
    · rw [Category.assoc, hκc m, ← Category.assoc, hκc n, Category.assoc, hmul]
    · rw [P.Base_comp, hκb m, hκb n, P.Base_id, Category.comp_id]
  · intro n
    exact ⟨hκb n, ⟨⟨prop_1_4_i P _ (fun Y _ => hiso Y), hdivκ n⟩,
      by show IsIso (P.Base (κ n)); rw [hκb n, P.Base_id]; infer_instance⟩⟩

/-- ★★★**Frobenius-trivial 対象の間では、底の任意の射が pull-back 射に持ち上がる**。

原文 (FrdI p.99):
> tion to D (respectively, by its Frobenius degree). Moreover, by assertion (iii), it fol-
-/
theorem exists_pullBack_frobTrivial (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A A' : C} (hA : IsFrobeniusTrivial P A) (hA' : IsFrobeniusTrivial P A')
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj A').base) :
    ∃ π : A ⟶ A', IsPullBack P π ∧ P.Base π = θ := by
  obtain ⟨A₀, π₀, ρ, hpb₀, hρ⟩ := exists_pullBack_over G A' θ
  have hA₀ := isFrobeniusTrivial_pullBack G hiso π₀ hpb₀ hA'
  haveI : IsIso ρ.inv := inferInstance
  obtain ⟨ι, hι⟩ := frobTrivial_iso_of_baseIso G hiso hA₀ hA ρ.inv
  refine ⟨ι.hom ≫ π₀, IsPullBack.comp P (isPullBack_of_isIso P ι.hom) hpb₀, ?_⟩
  rw [P.Base_comp, hι, ← hρ, ← Category.assoc, Iso.inv_hom_id, Category.id_comp]

/-! ## ★3. 底変換 —— 類は `Φ(θ)` で写る -/

/-- ★★★★**[FrdI] Theorem 5.1, (ii) の底変換** ——
`θ : A_𝒟 → A'_𝒟` に沿って `Pic` の対象を引き戻すと、類はちょうど `Φ(θ)` で写る。 -/
theorem picPullBack (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A A' : C} (hA : IsFrobeniusTrivial P A) (hA' : IsFrobeniusTrivial P A')
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj A').base)
    (Z' : PicObj P A') {c' : Gp (Φ.val (P.toElem.obj A').base)} (hc' : Z'.HasCls c') :
    ∃ (Z'' : PicObj P A) (p : Z''.obj ⟶ Z'.obj), IsPullBack P p ∧
      P.Base p = Z''.iso.hom ≫ θ ≫ Z'.iso.inv ∧
      Z''.HasCls (gpMap _ (Φ.map θ) c') := by
  obtain ⟨π, hπpb, hπb⟩ := exists_pullBack_frobTrivial G hiso hA hA' θ
  obtain ⟨X', φ', ψ', hsφ', hsψ', hbsp, hcls'⟩ := hc'
  haveI hbφ' : IsIso (P.Base φ') := hsφ'.2
  obtain ⟨B'', p, ρ'', hppb, hρ''⟩ := exists_pullBack_over G Z'.obj (θ ≫ Z'.iso.inv)
  obtain ⟨X'', q, ρX, hqpb, hρX⟩ := exists_pullBack_over G X' (θ ≫ inv (P.Base φ'))
  haveI hbρX : IsIso ρX.hom := ρX.isIso_hom
  obtain ⟨φ'', hφ''c, hφ''b⟩ := pullBack_lift hπpb (q ≫ φ') ρX.hom (by
    rw [P.Base_comp, ← hρX, hπb, Category.assoc, Category.assoc, IsIso.inv_hom_id,
      Category.comp_id])
  obtain ⟨ψ'', hψ''c, hψ''b⟩ := pullBack_lift hppb (q ≫ ψ') (ρX.hom ≫ ρ''.inv) (by
    rw [P.Base_comp, ← hρX, ← hρ'', Category.assoc, Category.assoc,
      picSpan_alpha (Z := Z') hsφ' hbsp, Category.assoc, ← Category.assoc ρ''.inv,
      Iso.inv_hom_id, Category.id_comp])
  have hdegπ : P.degFr π = 1 := (G.core.pullBackLB π hπpb).2
  have hdegp : P.degFr p = 1 := (G.core.pullBackLB p hppb).2
  have hdegq : P.degFr q = 1 := (G.core.pullBackLB q hqpb).2
  have hdivπ : P.Div π = 0 := (G.core.pullBackLB π hπpb).1.2
  have hdivp : P.Div p = 0 := (G.core.pullBackLB p hppb).1.2
  have hdivq : P.Div q = 0 := (G.core.pullBackLB q hqpb).1.2
  have hdegφ'' : P.degFr φ'' = 1 := by
    have h := congrArg P.degFr hφ''c
    rw [P.degFr_comp, P.degFr_comp, hdegπ, hdegq, hsφ'.1, one_mul, one_mul] at h
    exact h
  have hdegψ'' : P.degFr ψ'' = 1 := by
    have h := congrArg P.degFr hψ''c
    rw [P.degFr_comp, P.degFr_comp, hdegp, hdegq, hsψ'.1, one_mul, one_mul] at h
    exact h
  have hsφ'' : IsPreStep P φ'' := ⟨hdegφ'', by
    show IsIso (P.Base φ''); rw [hφ''b]; exact ρX.isIso_hom⟩
  have hsψ'' : IsPreStep P ψ'' := ⟨hdegψ'', by
    show IsIso (P.Base ψ''); rw [hψ''b]; exact (ρX ≪≫ ρ''.symm).isIso_hom⟩
  haveI hbφ'' : IsIso (P.Base φ'') := hsφ''.2
  have hdφ'' : P.Div φ'' = Φ.map (P.Base q) (P.Div φ') := by
    have h := congrArg P.Div hφ''c
    rw [P.Div_comp, P.Div_comp, hdivπ, hdivq, map_zero, smul_zero, add_zero, zero_add,
      hdegπ] at h
    rw [← h]; show _ = ((1 : ℕ+) : ℕ) • P.Div φ''
    rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  have hdψ'' : P.Div ψ'' = Φ.map (P.Base q) (P.Div ψ') := by
    have h := congrArg P.Div hψ''c
    rw [P.Div_comp, P.Div_comp, hdivp, hdivq, map_zero, smul_zero, add_zero, zero_add,
      hdegp] at h
    rw [← h]; show _ = ((1 : ℕ+) : ℕ) • P.Div ψ''
    rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  have hqb : P.Base q ≫ P.Base φ' = ρX.hom ≫ θ := by
    rw [← hρX, Category.assoc, Category.assoc, IsIso.inv_hom_id, Category.comp_id]
  have hinvb : @inv _ _ _ _ (P.Base φ'') hbφ'' = inv ρX.hom := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [hφ''b, IsIso.hom_inv_id]
  refine ⟨⟨B'', ρ''⟩, p, hppb, hρ''.symm, X'', φ'', ψ'', hsφ'', hsψ'', ?_, ?_⟩
  · show P.Base ψ'' ≫ ρ''.hom = P.Base φ''
    rw [hψ''b, hφ''b, Category.assoc, Iso.inv_hom_id, Category.comp_id]
  · rw [spanCls_eq, hinvb, hdφ'', hdψ'', ← gpMap_toGp, ← gpMap_toGp, ← map_sub,
      ← hcls', ← gpMap_base_spanCls φ' hsφ'.2 ψ',
      gpMap_phi_comp Φ (P.Base q) (P.Base φ'), hqb, ← gpMap_phi_comp Φ ρX.hom θ]
    exact gpMap_phi_inv_right ρX.hom _

/-! ## ★4. pre-step で類は `Div` ぶんずれる -/

/-- ★★**pre-step で `Pic` の類は `Div` だけずれる**。 -/
theorem PicObj.hasCls_preStep {A : C} {Z : PicObj P A}
    {c : Gp (Φ.val (P.toElem.obj A).base)} (h : Z.HasCls c)
    {V : C} (β₀ : Z.obj ⟶ V) (hs : IsPreStep P β₀) :
    (⟨V, (@asIso _ _ _ _ (P.Base β₀) hs.2).symm ≪≫ Z.iso⟩ : PicObj P A).HasCls
      (c + toGp _ (Φ.map Z.iso.inv (P.Div β₀))) := by
  haveI hbβ : IsIso (P.Base β₀) := hs.2
  obtain ⟨X, φ, ψ, hsφ, hsψ, hb, rfl⟩ := h
  haveI hbφ : IsIso (P.Base φ) := hsφ.2
  refine ⟨X, φ, ψ ≫ β₀, hsφ, IsPreStep.comp P hsψ hs, ?_, ?_⟩
  · show P.Base (ψ ≫ β₀) ≫ (inv (P.Base β₀) ≫ Z.iso.hom) = P.Base φ
    rw [P.Base_comp, Category.assoc, ← Category.assoc (P.Base β₀), IsIso.hom_inv_id,
      Category.id_comp, hb]
  · rw [spanCls_eq, spanCls_eq, Div_comp_preStep _ _ hs.1, toGp_add, ← gpMap_toGp,
      show (gpMap (Φ.val (P.toElem.obj Z.obj).base) (Φ.map (P.Base ψ))
            (toGp (Φ.val (P.toElem.obj Z.obj).base) (P.Div β₀))
          + toGp (Φ.val (P.toElem.obj X).base) (P.Div ψ))
          - toGp (Φ.val (P.toElem.obj X).base) (P.Div φ)
        = (toGp (Φ.val (P.toElem.obj X).base) (P.Div ψ)
            - toGp (Φ.val (P.toElem.obj X).base) (P.Div φ))
          + gpMap (Φ.val (P.toElem.obj Z.obj).base) (Φ.map (P.Base ψ))
            (toGp (Φ.val (P.toElem.obj Z.obj).base) (P.Div β₀)) from by abel,
      map_add, gpMap_phi_comp, picSpan_alpha (Z := Z) hsφ hb, gpMap_toGp]

/-! ## ★5. 同底の場合 -/

/-- ★★★★**[FrdI] Theorem 5.1, (ii) の同底版**(類の関係式 ⟹ 射の存在)。 -/
theorem picHom_sameBase_of_mem (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A : C} (hA : IsFrobeniusTrivial P A)
    (Z W : PicObj P A) {c cw : Gp (Φ.val (P.toElem.obj A).base)}
    (hc : Z.HasCls c) (hcw : W.HasCls cw)
    (d : ℕ+) (z : Φ.val (P.toElem.obj Z.obj).base)
    (heq : (d : ℕ) • c + toGp _ (Φ.map Z.iso.inv z) - cw ∈ phiBiratAt P G A) :
    ∃ φ : Z.obj ⟶ W.obj, P.degFr φ = d ∧
      P.Base φ = Z.iso.hom ≫ W.iso.inv ∧ P.Div φ = z := by
  obtain ⟨U, γ, hγ, hγdeg⟩ := G.core.frobDegSurj Z.obj d
  haveI hbγ : IsIso (P.Base γ) := hγ.2
  have hcU := PicObj.hasCls_frobType G hA hc γ hγ
  rw [hγdeg] at hcU
  obtain ⟨V, β₀, hβc, hβs, hβdiv⟩ := exists_coaPreStep_div G U (Φ.map (inv (P.Base γ)) z)
  haveI hbβ : IsIso (P.Base β₀) := hβs.2
  have hcV := PicObj.hasCls_preStep hcU β₀ hβs
  have hkey : Φ.map (Z.iso.inv ≫ P.Base γ) (P.Div β₀) = Φ.map Z.iso.inv z := by
    rw [hβdiv, Φ.map_comp, phi_map_comp_inv]
  rw [show ((⟨U, (@asIso _ _ _ _ (P.Base γ) hγ.2).symm ≪≫ Z.iso⟩ : PicObj P A).iso.inv)
      = Z.iso.inv ≫ P.Base γ from rfl, hkey] at hcV
  obtain ⟨ι, hι⟩ := PicObj.rel_of_cls_sub_mem G hiso hcV hcw heq
  refine ⟨γ ≫ β₀ ≫ ι.hom, ?_, ?_, ?_⟩
  · rw [P.degFr_comp, P.degFr_comp, degFr_of_isIso P ι.hom, hβs.1, hγdeg, one_mul, one_mul]
  · have hιb : P.Base ι.hom ≫ W.iso.hom
        = inv (P.Base β₀) ≫ inv (P.Base γ) ≫ Z.iso.hom := hι
    have h2 : P.Base ι.hom = (inv (P.Base β₀) ≫ inv (P.Base γ) ≫ Z.iso.hom) ≫ W.iso.inv := by
      rw [← hιb, Category.assoc, Iso.hom_inv_id, Category.comp_id]
    rw [P.Base_comp, P.Base_comp, h2]
    simp
  · rw [P.Div_comp, P.Div_comp, show P.Div ι.hom = 0 from isIsometric_of_isIso P ι.hom,
      map_zero, zero_add, degFr_of_isIso P ι.hom,
      show P.Div γ = 0 from hγ.1.2, smul_zero, add_zero]
    rw [show (((1 : ℕ+) : ℕ) • P.Div β₀) = P.Div β₀ from by
      rw [show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul], hβdiv, phi_map_comp_inv]

/-- ★★★★**[FrdI] Theorem 5.1, (ii) の同底版**(射の存在 ⟹ 類の関係式)。 -/
theorem picHom_sameBase_mem_of (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A : C} (hA : IsFrobeniusTrivial P A)
    (Z W : PicObj P A) {c cw : Gp (Φ.val (P.toElem.obj A).base)}
    (hc : Z.HasCls c) (hcw : W.HasCls cw)
    (d : ℕ+) (z : Φ.val (P.toElem.obj Z.obj).base)
    (φ : Z.obj ⟶ W.obj) (hdeg : P.degFr φ = d)
    (hbase : P.Base φ = Z.iso.hom ≫ W.iso.inv) (hdiv : P.Div φ = z) :
    (d : ℕ) • c + toGp _ (Φ.map Z.iso.inv z) - cw ∈ phiBiratAt P G A := by
  obtain ⟨U, V, γ, β₀, α, hfac, hγ, hβs, hαpb⟩ := G.core.arbFactor φ
  haveI hbγ : IsIso (P.Base γ) := hγ.2
  haveI hbβ : IsIso (P.Base β₀) := hβs.2
  have hdegα : P.degFr α = 1 := (G.core.pullBackLB α hαpb).2
  have hdivα : P.Div α = 0 := (G.core.pullBackLB α hαpb).1.2
  have hdivγ : P.Div γ = 0 := hγ.1.2
  have hdegγ : P.degFr γ = d := by
    have h := congrArg P.degFr hfac
    rw [hdeg, P.degFr_comp, P.degFr_comp, hdegα, hβs.1, one_mul, one_mul] at h
    exact h.symm
  haveI hbφ : IsIso (P.Base φ) := by rw [hbase]; infer_instance
  haveI hbα : IsIso (P.Base α) := by
    have h := congrArg P.Base hfac
    rw [P.Base_comp, P.Base_comp] at h
    haveI h2 : IsIso (P.Base γ ≫ P.Base β₀ ≫ P.Base α) := by rw [← h]; exact hbφ
    haveI h3 : IsIso (P.Base β₀ ≫ P.Base α) :=
      IsIso.of_isIso_comp_left (P.Base γ) (P.Base β₀ ≫ P.Base α)
    exact IsIso.of_isIso_comp_left (P.Base β₀) (P.Base α)
  haveI hαiso : IsIso α := isIso_of_pullBack_base_isIso hαpb hbα
  have hdivβ : P.Div β₀ = Φ.map (inv (P.Base γ)) z := by
    have h := congrArg P.Div hfac
    rw [hdiv, P.Div_comp, P.Div_comp, hdivα, hdivγ, map_zero, zero_add, smul_zero, add_zero,
      hdegα] at h
    rw [h, show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul, phi_map_inv_comp]
  have hcU := PicObj.hasCls_frobType G hA hc γ hγ
  rw [hdegγ] at hcU
  have hcV := PicObj.hasCls_preStep hcU β₀ hβs
  have hkey : Φ.map (Z.iso.inv ≫ P.Base γ) (P.Div β₀) = Φ.map Z.iso.inv z := by
    rw [hdivβ, Φ.map_comp, phi_map_comp_inv]
  rw [show ((⟨U, (@asIso _ _ _ _ (P.Base γ) hγ.2).symm ≪≫ Z.iso⟩ : PicObj P A).iso.inv)
      = Z.iso.inv ≫ P.Base γ from rfl, hkey] at hcV
  have hrel : PicObj.Rel
      (⟨V, (@asIso _ _ _ _ (P.Base β₀) hβs.2).symm ≪≫
        ((@asIso _ _ _ _ (P.Base γ) hγ.2).symm ≪≫ Z.iso)⟩ : PicObj P A) W := by
    refine ⟨asIso α, ?_⟩
    show P.Base α ≫ W.iso.hom = inv (P.Base β₀) ≫ inv (P.Base γ) ≫ Z.iso.hom
    have h := congrArg P.Base hfac
    rw [P.Base_comp, P.Base_comp, hbase] at h
    have hα : P.Base α = inv (P.Base β₀) ≫ inv (P.Base γ) ≫ Z.iso.hom ≫ W.iso.inv := by
      rw [h]; simp
    rw [hα]; simp
  exact PicObj.cls_sub_mem G hiso (PicObj.hasCls_of_rel hrel hcV) hcw

/-! ## ★6. ★★★★`Theorem 5.1, (ii)` -/

/-- ★★★★**[FrdI] Theorem 5.1, (ii)** —— 射の存在と類の関係式は同値。

原文 (FrdI p.97):
> [where, by abuse of notation, we denote by z|AD the image of
-/
theorem thm_5_1_ii (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A A' : C} (hA : IsFrobeniusTrivial P A) (hA' : IsFrobeniusTrivial P A')
    (Z : PicObj P A) (Z' : PicObj P A')
    {c : Gp (Φ.val (P.toElem.obj A).base)} {c' : Gp (Φ.val (P.toElem.obj A').base)}
    (hc : Z.HasCls c) (hc' : Z'.HasCls c')
    (θ : (P.toElem.obj A).base ⟶ (P.toElem.obj A').base) (d : ℕ+)
    (z : Φ.val (P.toElem.obj Z.obj).base) :
    (∃ φ : Z.obj ⟶ Z'.obj, P.degFr φ = d ∧
        P.Base φ = Z.iso.hom ≫ θ ≫ Z'.iso.inv ∧ P.Div φ = z)
      ↔ (d : ℕ) • c + toGp _ (Φ.map Z.iso.inv z) - gpMap _ (Φ.map θ) c'
          ∈ phiBiratAt P G A := by
  obtain ⟨Z'', p, hppb, hpb, hc''⟩ := picPullBack G hiso hA hA' θ Z' hc'
  have hdegp : P.degFr p = 1 := (G.core.pullBackLB p hppb).2
  have hdivp : P.Div p = 0 := (G.core.pullBackLB p hppb).1.2
  constructor
  · rintro ⟨φ, hdeg, hbase, hdivφ⟩
    obtain ⟨φ'', hφ''c, hφ''b⟩ := pullBack_lift hppb φ (Z.iso.hom ≫ Z''.iso.inv) (by
      rw [hbase, hpb, Category.assoc, ← Category.assoc Z''.iso.inv, Iso.inv_hom_id,
        Category.id_comp])
    refine picHom_sameBase_mem_of G hiso hA Z Z'' hc hc'' d z φ'' ?_ hφ''b ?_
    · have h := congrArg P.degFr hφ''c
      rw [P.degFr_comp, hdegp, one_mul, hdeg] at h
      exact h
    · have h := congrArg P.Div hφ''c
      rw [P.Div_comp, hdivp, map_zero, zero_add, hdegp, hdivφ] at h
      rw [← h, show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul]
  · intro hmem
    obtain ⟨φ'', hdeg'', hbase'', hdiv''⟩ :=
      picHom_sameBase_of_mem G hiso hA Z Z'' hc hc'' d z hmem
    refine ⟨φ'' ≫ p, ?_, ?_, ?_⟩
    · rw [P.degFr_comp, hdegp, one_mul, hdeg'']
    · rw [P.Base_comp, hbase'', hpb, Category.assoc, ← Category.assoc Z''.iso.inv,
        Iso.inv_hom_id, Category.id_comp]
    · rw [P.Div_comp, hdivp, map_zero, zero_add, hdegp,
        show ((1 : ℕ+) : ℕ) = 1 from rfl, one_smul, hdiv'']

/-- ★★**[FrdI] Theorem 5.1, (ii) の「Moreover」** —— 底・零因子・次数が一致すれば
`𝔽_Φ`(= `𝒞^un-tr`)の像は等しい —— すなわち**単元同値類は一意**。

原文 (FrdI p.97):
> Moreover, if such a morphism exists, then its unit-equivalence class [i.e., its image
-/
theorem thm_5_1_ii_uniq {A B : C} (φ φ' : A ⟶ B)
    (hb : P.Base φ = P.Base φ') (hz : P.Div φ = P.Div φ') (hd : P.degFr φ = P.degFr φ') :
    P.toElem.map φ = P.toElem.map φ' :=
  ElemFrobCat.Hom.ext hb hz hd

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.1, (ii)`。 -/
def thm_5_1_ii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 96, item := "Theorem 5.1, (ii)",
    sectionId := "frdi-thm-5-1" }

end ABC3.Found.FrdI
