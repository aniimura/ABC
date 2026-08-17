import ABC3.Found.FrdI.Def45
import ABC3.Found.FrdI.PlBkShuffle
import ABC3.Found.FrdI.Remark311

/-!
# [FrdI] Example 4.3 —— 右手と左手の同型は一致しない

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.81–82。

★★**原文の構成**:

- 底圏 `𝒟` は**射が 1 本しかない圏**。
- `Φ` は `𝒟` の唯一の対象での値が `ℚ≥0` の単系。
- `𝒞` の対象は `ℚ` の元、射 `a → b` は **`d · a ≤ b` なる `d ∈ ℕ≥1`**。
  合成は `ℕ≥1` の掛け算。
- 関手 `𝒞 → 𝔽_Φ` は `φ : a → b` に零因子 `b − deg_Fr(φ)·a ∈ ℚ≥0` と
  Frobenius 次数 `deg_Fr(φ)` を対応させる。

★★★**この例が具体的である理由**: 射は**次数 `d` だけで決まる**
(`d · a ≤ b` は命題なので、`Hom(a,b)` は `ℕ≥1` の部分集合)。
したがって `Definition 1.3` の 21 条が**すべて `ℚ` の不等式の計算**になる。

原文 (FrdI p.81):
> Independence of Right-hand and Left-hand Isomor-
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe u v

/-! ## ★1. `ℚ≥0` は divisorial

★`sharp` は「和が `0` なら両方 `0`」、`integral` は簡約性、
`saturated` は `n·b ≤ n·a → b ≤ a`、`of characteristic type` は sharp から。 -/

/-- ★`ℚ≥0` は **sharp** —— 可逆元は `0` だけ。 -/
theorem isSharp_nnrat : IsSharp ℚ≥0 := by
  intro a ha
  obtain ⟨u, rfl⟩ := ha
  have h := u.val_neg
  exact (add_eq_zero.mp h).1

/-- ★`ℚ≥0` は **integral** —— 簡約的だから。 -/
theorem isIntegralMonoid_nnrat : IsIntegralMonoid ℚ≥0 :=
  isIntegralMonoid_of_isCancelAdd ℚ≥0

/-- ★`ℚ≥0` の `≼` は順序そのもの。 -/
theorem mle_nnrat_iff {a b : ℚ≥0} : MLe a b ↔ a ≤ b := by
  constructor
  · rintro ⟨c, rfl⟩
    exact le_self_add
  · intro h
    exact le_iff_exists_add.mp h |>.imp fun c hc => hc.symm

/-- ★`ℚ≥0` は **saturated**。 -/
theorem isSaturatedMonoid_nnrat : IsSaturatedMonoid ℚ≥0 := by
  rw [isSaturatedMonoid_iff_mle ℚ≥0 inferInstance]
  intro a b n hn hle
  rw [mle_nnrat_iff] at hle ⊢
  have hn' : (0 : ℚ≥0) < (n : ℚ≥0) := by exact_mod_cast hn
  rw [nsmul_eq_mul, nsmul_eq_mul] at hle
  exact le_of_mul_le_mul_left hle hn'

/-- ★`ℚ≥0` は **pre-divisorial**。 -/
theorem isPreDivisorial_nnrat : IsPreDivisorial ℚ≥0 :=
  ⟨isIntegralMonoid_nnrat, isSaturatedMonoid_nnrat,
    isOfCharacteristicType_of_isSharp ℚ≥0 isSharp_nnrat⟩

/-- ★★`ℚ≥0` は **divisorial**。 -/
theorem isDivisorial_nnrat : IsDivisorial ℚ≥0 :=
  ⟨isPreDivisorial_nnrat, isSharp_nnrat⟩

/-! ## ★2. 底圏と `Φ` —— 一射圏の上の定数 `ℚ≥0`

原文 (FrdI p.82):
> Let D be a one-morphism category; Φ the monoid on D whose value on the unique
-/

/-- ★★`Example 4.3` の `Φ` —— 一射圏の上の定数 `ℚ≥0`。 -/
abbrev Phi43 : MonoidOn.{0, 0, 0} (Discrete PUnit) := constPhi ℚ≥0

/-! ## ★3. `𝒞` —— 対象は `ℚ`、射 `a → b` は `d·a ≤ b` なる `d ∈ ℕ≥1`

原文 (FrdI p.82):
> morphisms is defined by multiplication of elements of N
-/

/-- ★★`Example 4.3` の `𝒞` の対象 —— `ℚ` の元(型シノニム)。 -/
def Ex43 : Type := ℚ

/-- ★対象の台となる有理数。 -/
def Ex43.val (a : Ex43) : ℚ := a

/-- ★有理数から対象へ。 -/
def Ex43.mk (q : ℚ) : Ex43 := q

@[simp] theorem Ex43.val_mk (q : ℚ) : (Ex43.mk q).val = q := rfl

@[ext] theorem Ex43.ext {a b : Ex43} (h : a.val = b.val) : a = b := h

/-- ★★`𝒞` の射 —— `d · a ≤ b` なる `d ∈ ℕ≥1`。 -/
@[ext]
structure Ex43Hom (a b : Ex43) where
  /-- Frobenius 次数 -/
  deg : ℕ+
  /-- `d · a ≤ b` -/
  cond : ((deg : ℕ) : ℚ) * a.val ≤ b.val

/-- ★★合成は `ℕ≥1` の掛け算。 -/
instance : Category Ex43 where
  Hom := Ex43Hom
  id a := ⟨1, by simp⟩
  comp {a b c} f g :=
    ⟨g.deg * f.deg, by
      have hg : (0 : ℚ) < ((g.deg : ℕ) : ℚ) := by
        exact_mod_cast g.deg.pos
      calc (((g.deg * f.deg : ℕ+) : ℕ) : ℚ) * a.val
          = ((g.deg : ℕ) : ℚ) * (((f.deg : ℕ) : ℚ) * a.val) := by push_cast; ring
        _ ≤ ((g.deg : ℕ) : ℚ) * b.val := by
            exact mul_le_mul_of_nonneg_left f.cond (le_of_lt hg)
        _ ≤ c.val := g.cond⟩
  id_comp f := Ex43Hom.ext (by show f.deg * 1 = f.deg; rw [mul_one])
  comp_id f := Ex43Hom.ext (by show 1 * f.deg = f.deg; rw [one_mul])
  assoc f g h := Ex43Hom.ext (by
    show h.deg * (g.deg * f.deg) = h.deg * g.deg * f.deg
    rw [mul_assoc])

@[simp] theorem Ex43.id_deg (a : Ex43) : (𝟙 a : Ex43Hom a a).deg = 1 := rfl

@[simp] theorem Ex43.comp_deg {a b c : Ex43} (f : a ⟶ b) (g : b ⟶ c) :
    (f ≫ g : Ex43Hom a c).deg = g.deg * f.deg := rfl

/-- ★★**射は次数だけで決まる**——`Hom` は `ℕ≥1` の部分集合。 -/
theorem Ex43.hom_ext {a b : Ex43} {f g : a ⟶ b} (h : f.deg = g.deg) : f = g :=
  Ex43Hom.ext h

/-! ## ★4. 関手 `𝒞 → 𝔽_Φ`

原文 (FrdI p.82):
> by assigning to a morphism φ : a
-/

/-- ★★零因子 `b − d·a ∈ ℚ≥0`。 -/
def ex43Div {a b : Ex43} (f : a ⟶ b) : ℚ≥0 :=
  ⟨b.val - ((f.deg : ℕ) : ℚ) * a.val, by linarith [Ex43Hom.cond f]⟩

@[simp] theorem ex43Div_val {a b : Ex43} (f : a ⟶ b) :
    (ex43Div f : ℚ) = b.val - ((f.deg : ℕ) : ℚ) * a.val := rfl

/-- ★★★**`𝒞 → 𝔽_Φ`**。 -/
noncomputable def ex43ToElem : Ex43 ⥤ ElemFrobCat Phi43 where
  obj _ := constObj ℚ≥0
  map {a b} f := ⟨𝟙 _, ex43Div f, f.deg⟩
  map_id a := by
    refine ElemFrobCat.Hom.ext rfl ?_ rfl
    refine NNRat.coe_injective ?_
    show ((ex43Div (𝟙 a) : ℚ≥0) : ℚ) = ((0 : ℚ≥0) : ℚ)
    rw [ex43Div_val, Ex43.id_deg]
    push_cast
    ring
  map_comp {a b c} f g := by
    refine ElemFrobCat.Hom.ext (Subsingleton.elim _ _) ?_ rfl
    rw [ElemFrobCat.div_comp]
    refine NNRat.coe_injective ?_
    rw [MonoidOn.const_map]
    dsimp only
    show ((ex43Div (f ≫ g) : ℚ≥0) : ℚ)
      = ((ex43Div g + ((g.deg : ℕ+) : ℕ) • ex43Div f : ℚ≥0) : ℚ)
    rw [nsmul_eq_mul]
    push_cast [ex43Div_val, Ex43.comp_deg]
    ring

/-! ## ★5. `𝒞` は連結かつ totally epimorphic —— pre-Frobenioid 構造

原文 (FrdI p.82):
> of C as the Frobenius degree of the morphism. Thus, we obtain a natural functor
-/

/-- ★★**`𝒞` は totally epimorphic** —— 射は次数だけで決まり、`ℕ≥1` は簡約的。 -/
theorem ex43_totEpi : IsTotallyEpimorphic Ex43 := by
  intro a b f
  refine ⟨fun {c} g h hgh => Ex43.hom_ext ?_⟩
  have hd : (f ≫ g : Ex43Hom a c).deg = (f ≫ h : Ex43Hom a c).deg :=
    congrArg (fun t : a ⟶ c => Ex43Hom.deg t) hgh
  rw [Ex43.comp_deg, Ex43.comp_deg] at hd
  exact mul_right_cancel hd

/-- ★底圏は totally epimorphic(射が 1 本しかない)。 -/
theorem d43_totEpi : IsTotallyEpimorphic (Discrete PUnit) :=
  fun _ _ _ => ⟨fun g h _ => Subsingleton.elim g h⟩

/-- ★底圏は連結。 -/
theorem d43_connected : IsConnected (Discrete PUnit) := by
  haveI : Nonempty (Discrete PUnit) := ⟨⟨⟨⟩⟩⟩
  refine zigzag_isConnected (fun X Y => ?_)
  obtain ⟨x⟩ := X
  obtain ⟨y⟩ := Y
  exact Relation.ReflTransGen.refl

/-- ★★**`𝒞` は連結** —— `a → max(a,b) ← b`。 -/
theorem ex43_connected : IsConnected Ex43 := by
  haveI : Nonempty Ex43 := ⟨Ex43.mk 0⟩
  refine zigzag_isConnected (fun X Y => ?_)
  have hX : X ⟶ Ex43.mk (max X.val Y.val) :=
    ⟨1, by simpa using le_max_left X.val Y.val⟩
  have hY : Y ⟶ Ex43.mk (max X.val Y.val) :=
    ⟨1, by simpa using le_max_right X.val Y.val⟩
  exact Relation.ReflTransGen.trans (Relation.ReflTransGen.single (Or.inl ⟨hX⟩))
    (Relation.ReflTransGen.single (Or.inr ⟨hY⟩))

/-- ★★★**`Example 4.3` の pre-Frobenioid 構造**。 -/
noncomputable def ex43P : PreFrobenioid Ex43 Phi43 where
  toElem := ex43ToElem
  divisorial _ := isDivisorial_nnrat
  totEpiC := ex43_totEpi
  totEpiD := d43_totEpi
  connectedC := ex43_connected
  connectedD := d43_connected

@[simp] theorem ex43P_degFr {a b : Ex43} (f : a ⟶ b) : ex43P.degFr f = f.deg := rfl

@[simp] theorem ex43P_Div {a b : Ex43} (f : a ⟶ b) : ex43P.Div f = ex43Div f := rfl

@[simp] theorem ex43P_Base {a b : Ex43} (f : a ⟶ b) :
    ex43P.Base f = 𝟙 (constObj ℚ≥0).base := rfl

/-! ## ★6. 射の型の判定

原文 (FrdI p.82):
> all morphisms of C are base-isomorphisms; all pull-back morphisms of C are iso-

★★底圏の射は 1 本しかないので **`Base` はつねに同型**、
したがって `pre-step ⟺ 次数 1`、`Frobenius 型 ⟺ 等長`。 -/

theorem ex43_isBaseIso {a b : Ex43} (f : a ⟶ b) : IsBaseIsomorphism ex43P f := by
  show IsIso (ex43P.Base f)
  exact ⟨𝟙 _, Subsingleton.elim _ _, Subsingleton.elim _ _⟩

theorem ex43_isLinear_iff {a b : Ex43} (f : a ⟶ b) : IsLinear ex43P f ↔ f.deg = 1 :=
  Iff.rfl

theorem ex43_isometric_iff {a b : Ex43} (f : a ⟶ b) :
    IsIsometric ex43P f ↔ b.val = ((f.deg : ℕ) : ℚ) * a.val := by
  constructor
  · intro h
    have h1 : ex43Div f = 0 := h
    have h2 : ((ex43Div f : ℚ≥0) : ℚ) = ((0 : ℚ≥0) : ℚ) := congrArg _ h1
    rw [ex43Div_val] at h2
    simp only [NNRat.coe_zero] at h2
    linarith
  · intro h
    show ex43Div f = 0
    refine NNRat.coe_injective ?_
    rw [ex43Div_val, h]
    simp

theorem ex43_isPreStep_iff {a b : Ex43} (f : a ⟶ b) : IsPreStep ex43P f ↔ f.deg = 1 :=
  ⟨fun h => h.1, fun h => ⟨h, ex43_isBaseIso f⟩⟩

/-- ★★★**等長 pre-step は同型** —— 等長は `b = d·a`、pre-step は `d = 1`、
合わせて `b = a` になり、射は次数だけで決まるので恒等射である。 -/
theorem ex43_isIso_of_isometric_preStep {a b : Ex43} (f : a ⟶ b)
    (hm : IsIsometric ex43P f) (hs : IsPreStep ex43P f) : IsIso f := by
  have hd : f.deg = 1 := hs.1
  have hb : b.val = a.val := by
    have h2 := (ex43_isometric_iff f).mp hm
    rw [hd] at h2
    simpa using h2
  have hab : a = b := Ex43.ext hb.symm
  subst hab
  have hfid : f = 𝟙 a := Ex43.hom_ext (by rw [hd, Ex43.id_deg])
  rw [hfid]
  infer_instance

/-- ★★**すべての対象が isotropic**。 -/
theorem ex43_isotropicType : IsOfIsotropicType ex43P :=
  fun _ _ φ hm hs => ex43_isIso_of_isometric_preStep φ hm hs

/-- ★★**すべての射が co-angular** —— 真ん中が等長 pre-step なら恒等射だから。 -/
theorem ex43_coAngular {a b : Ex43} (f : a ⟶ b) : IsCoAngular ex43P f :=
  fun _ _ _ β _ _ _ hm hs _ => ex43_isIso_of_isometric_preStep β hm hs

theorem ex43_isFrobType_iff {a b : Ex43} (f : a ⟶ b) :
    IsFrobeniusType ex43P f ↔ IsIsometric ex43P f :=
  ⟨fun h => h.1.2, fun h => ⟨⟨ex43_coAngular f, h⟩, ex43_isBaseIso f⟩⟩

theorem ex43_isLBInv_iff {a b : Ex43} (f : a ⟶ b) :
    IsLBInvertible ex43P f ↔ IsIsometric ex43P f :=
  ⟨fun h => h.2, fun h => ⟨ex43_coAngular f, h⟩⟩

/-- ★`ℕ≥1` で `d · e = 1` なら `d = 1`。 -/
theorem ex43_pnat_eq_one {d e : ℕ+} (h : d * e = 1) : d = 1 := by
  have hcoe : ((d : ℕ+) : ℕ) * ((e : ℕ+) : ℕ) = 1 := by rw [← PNat.mul_coe, h]; rfl
  exact PNat.coe_injective (Nat.dvd_one.mp ⟨_, hcoe.symm⟩)

/-- ★★★**pull-back 射は同型** —— `b` から `𝟙 b` を持ち上げると
`f.deg · h.deg = 1` かつ `h.deg · b ≤ a` が出る。 -/
theorem ex43_isIso_of_isPullBack {a b : Ex43} (f : a ⟶ b) (hf : IsPullBack ex43P f) :
    IsIso f := by
  obtain ⟨h, ⟨hh1, -⟩, -⟩ :=
    IsPullBack.lift ex43P hf b (𝟙 b) (𝟙 _) (Subsingleton.elim _ _)
  have hdeg : f.deg * h.deg = 1 := by
    have h2 := congrArg (fun t : b ⟶ b => Ex43Hom.deg t) hh1
    rw [Ex43.comp_deg, Ex43.id_deg] at h2
    exact h2
  have hf1 : f.deg = 1 := ex43_pnat_eq_one hdeg
  have hh1' : h.deg = 1 := ex43_pnat_eq_one (by rw [mul_comm]; exact hdeg)
  have hab : a = b := by
    refine Ex43.ext (le_antisymm ?_ ?_)
    · have := Ex43Hom.cond f
      rw [hf1] at this
      simpa using this
    · have := Ex43Hom.cond h
      rw [hh1'] at this
      simpa using this
  subst hab
  have hfid : f = 𝟙 a := Ex43.hom_ext (by rw [hf1, Ex43.id_deg])
  rw [hfid]
  infer_instance

/-! ## ★7. 対象の等式に落ちる 2 つの補題

★★`Example 4.3` では**同型は恒等射しかない**(圏は skeletal)。
その源は「等長 pre-step も pull-back も対象を動かさない」ことである。 -/

instance : Subsingleton (Discrete PUnit) :=
  ⟨fun a b => by obtain ⟨x⟩ := a; obtain ⟨y⟩ := b; rfl⟩

theorem ex43_eq_of_isometric_preStep {a b : Ex43} (f : a ⟶ b)
    (hm : IsIsometric ex43P f) (hs : IsPreStep ex43P f) : a = b := by
  have hd : f.deg = 1 := hs.1
  have h2 := (ex43_isometric_iff f).mp hm
  rw [hd] at h2
  have h3 : a.val = b.val := by simpa using h2.symm
  exact Ex43.ext h3

theorem ex43_eq_of_isPullBack {a b : Ex43} (f : a ⟶ b) (hf : IsPullBack ex43P f) :
    a = b := by
  haveI := ex43_isIso_of_isPullBack f hf
  exact ex43_eq_of_isometric_preStep f (isIsometric_of_isIso ex43P f)
    (isPreStep_of_isIso ex43P f)

theorem ex43_deg_of_isPullBack {a b : Ex43} (f : a ⟶ b) (hf : IsPullBack ex43P f) :
    f.deg = 1 := by
  haveI := ex43_isIso_of_isPullBack f hf
  exact isLinear_of_isIso ex43P f

/-- ★`𝒪^▷(a)` の元は恒等射だけ。 -/
theorem ex43_oTri_eq_id {a : Ex43} (α : End a) (hα : α ∈ OTri ex43P a) :
    α = 𝟙 a := Ex43.hom_ext (by rw [Ex43.id_deg]; exact hα.2)

/-! ## ★8. (i)(c) —— `(𝒞^pl-bk)_A → 𝒟_{A_𝒟}` は圏同値

★★pull-back 射は同型しかなく、同型は恒等射しかないので、
**両辺とも「対象 1 個・射 1 本」の圏**である。 -/

theorem ex43_plBkOver_faithful (A : Ex43) : (plBkOverFunctor ex43P A).Faithful where
  map_injective {Z W} {f g} _ := by
    refine Over.OverMorphism.ext (InducedWideCategory.Hom.ext ?_)
    exact Ex43.hom_ext (by
      rw [ex43_deg_of_isPullBack _ f.left.property,
        ex43_deg_of_isPullBack _ g.left.property])

theorem ex43_plBkOver_full (A : Ex43) : (plBkOverFunctor ex43P A).Full := by
  constructor
  intro Z W h
  have hZ : Z.left.obj = A := ex43_eq_of_isPullBack Z.hom.hom Z.hom.property
  have hW : W.left.obj = A := ex43_eq_of_isPullBack W.hom.hom W.hom.property
  have hval : Z.left.obj.val = W.left.obj.val := by rw [hZ, hW]
  have hle : (((1 : ℕ+) : ℕ) : ℚ) * Z.left.obj.val ≤ W.left.obj.val := by
    simpa using le_of_eq hval
  obtain ⟨g, hgdeg⟩ : ∃ g : Z.left.obj ⟶ W.left.obj, g.deg = 1 :=
    ⟨⟨1, hle⟩, rfl⟩
  haveI : IsIso g :=
    ex43_isIso_of_isometric_preStep g
      ((ex43_isometric_iff g).mpr (by rw [hgdeg]; simpa using hval.symm))
      ((ex43_isPreStep_iff g).mpr hgdeg)
  have hgc : g ≫ W.hom.hom = Z.hom.hom :=
    Ex43.hom_ext (by
      rw [Ex43.comp_deg, hgdeg, mul_one,
        ex43_deg_of_isPullBack _ W.hom.property,
        ex43_deg_of_isPullBack _ Z.hom.property])
  refine ⟨Over.homMk (show Z.left ⟶ W.left from ⟨g, isPullBack_of_isIso ex43P g⟩)
    (WideSubcategory.hom_ext _ hgc), ?_⟩
  exact Over.OverMorphism.ext (Subsingleton.elim _ _)

theorem ex43_plBkOver_essSurj (A : Ex43) : (plBkOverFunctor ex43P A).EssSurj := by
  refine ⟨fun T => ?_⟩
  refine ⟨Over.mk (show (⟨A⟩ : PlBk ex43P) ⟶ ⟨A⟩ from
    ⟨𝟙 A, isPullBack_of_isIso ex43P (𝟙 A)⟩), ⟨?_⟩⟩
  exact Over.isoMk (eqToIso (Subsingleton.elim _ _)) (Subsingleton.elim _ _)

theorem ex43_plBkEquiv (A : Ex43) : (plBkOverFunctor ex43P A).IsEquivalence :=
  ⟨ex43_plBkOver_faithful A, ex43_plBkOver_full A, ex43_plBkOver_essSurj A⟩

/-! ## ★9. `Definition 1.3` の 21 条

原文 (FrdI p.82):
> verifies immediately that C is a Frobenioid of isotropic type. Since D is clearly of
-/

/-- ★`0 ∈ ℚ` は Frobenius-trivial。 -/
theorem ex43_frobTrivial_zero : IsFrobeniusTrivial ex43P (Ex43.mk 0) := by
  refine ⟨{ toFun := fun n => (⟨n, by simp⟩ : End (Ex43.mk 0))
            map_one' := Ex43.hom_ext rfl
            map_mul' := fun x y => Ex43.hom_ext rfl }, fun n => rfl, fun n => ⟨?_, ?_⟩⟩
  · exact Subsingleton.elim _ _
  · refine (ex43_isFrobType_iff _).mpr ((ex43_isometric_iff _).mpr ?_)
    simp

/-- ★★★**`Example 4.3` の `𝒞` は `Definition 1.3` の 21 条を満たす**。 -/
theorem ex43_core : FrobenioidCore ex43P where
  baseSurj Y := ⟨Ex43.mk 0, ex43_frobTrivial_zero, ⟨eqToIso (Subsingleton.elim _ _)⟩⟩
  preStepSpan A B α _ :=
    ⟨Ex43.mk (min A.val B.val),
      ⟨1, by simpa using min_le_left A.val B.val⟩,
      ⟨1, by simpa using min_le_right A.val B.val⟩,
      (ex43_isPreStep_iff _).mpr rfl, (ex43_isPreStep_iff _).mpr rfl,
      Subsingleton.elim _ _⟩
  plBkEquiv A := ex43_plBkEquiv A
  frobDegSurj A n :=
    ⟨Ex43.mk (((n : ℕ+) : ℕ) * A.val), ⟨n, by simp⟩,
      (ex43_isFrobType_iff _).mpr ((ex43_isometric_iff _).mpr (by simp)), rfl⟩
  frobDegUniq A B E φ ψ hφ hψ hdeg := by
    have hB := (ex43_isometric_iff φ).mp ((ex43_isFrobType_iff φ).mp hφ)
    have hE := (ex43_isometric_iff ψ).mp ((ex43_isFrobType_iff ψ).mp hψ)
    have hdeg' : φ.deg = ψ.deg := hdeg
    have hBE : B = E := Ex43.ext (by rw [hB, hE, hdeg'])
    subst hBE
    exact ⟨𝟙 B, inferInstance, Ex43.hom_ext (by
      rw [Ex43.comp_deg, Ex43.id_deg, one_mul]; exact hdeg')⟩
  coAngularComp _ _ _ _ := ex43_coAngular _
  coAngularOfPreStep _ _ _ _ := ex43_coAngular _
  otriFwd φ _ _ α hα := by
    refine ⟨𝟙 _, ⟨(OTri ex43P _).one_mem, ?_⟩, ?_⟩
    · rw [ex43_oTri_eq_id α hα, Category.comp_id, Category.id_comp]
    · rintro y ⟨hy, -⟩
      exact ex43_oTri_eq_id y hy
  otriBwd φ _ _ β hβ := by
    refine ⟨𝟙 _, ⟨(OTri ex43P _).one_mem, ?_⟩, ?_⟩
    · rw [ex43_oTri_eq_id β hβ, Category.comp_id, Category.id_comp]
    · rintro y ⟨hy, -⟩
      exact ex43_oTri_eq_id y hy
  otriBase φ φ' _ _ _ _ _ α hα β hβ _ := by
    rw [ex43_oTri_eq_id α hα, ex43_oTri_eq_id β hβ, Category.comp_id, Category.id_comp]
  arbFactor {A B} φ :=
    ⟨Ex43.mk (((φ.deg : ℕ) : ℚ) * A.val), B,
      ⟨φ.deg, by simp⟩, ⟨1, by simpa using Ex43Hom.cond φ⟩, 𝟙 B,
      Ex43.hom_ext (by simp),
      (ex43_isFrobType_iff _).mpr ((ex43_isometric_iff _).mpr (by simp)),
      (ex43_isPreStep_iff _).mpr rfl,
      isPullBack_of_isIso ex43P (𝟙 B)⟩
  arbFactorUniq {A B} X Y X' Y' γ β α γ' β' α' heq hγ hβ hα hγ' hβ' hα' := by
    have hαd : α.deg = 1 := ex43_deg_of_isPullBack α hα
    have hα'd : α'.deg = 1 := ex43_deg_of_isPullBack α' hα'
    have hβd : β.deg = 1 := hβ.1
    have hβ'd : β'.deg = 1 := hβ'.1
    have hdeg : γ.deg = γ'.deg := by
      have hd := congrArg (fun t : A ⟶ B => Ex43Hom.deg t) heq
      rw [Ex43.comp_deg, Ex43.comp_deg, Ex43.comp_deg, Ex43.comp_deg,
        hαd, hα'd, hβd, hβ'd] at hd
      simpa using hd
    have hX : X = X' := by
      have h1 := (ex43_isometric_iff γ).mp ((ex43_isFrobType_iff γ).mp hγ)
      have h2 := (ex43_isometric_iff γ').mp ((ex43_isFrobType_iff γ').mp hγ')
      exact Ex43.ext (by rw [h1, h2, hdeg])
    have hY : Y = B := ex43_eq_of_isPullBack α hα
    have hY' : Y' = B := ex43_eq_of_isPullBack α' hα'
    subst hX
    subst hY
    subst hY'
    refine ⟨Iso.refl _, Iso.refl _, Ex43.hom_ext ?_, Ex43.hom_ext ?_, Ex43.hom_ext ?_⟩
    · rw [Ex43.comp_deg, hα'd, hαd]
      simp
    · rw [Ex43.comp_deg, Ex43.comp_deg, hβ'd, hβd]
      simp
    · rw [Ex43.comp_deg]
      simp [← hdeg]
  pullBackLB α hα := by
    haveI := ex43_isIso_of_isPullBack α hα
    exact ⟨⟨ex43_coAngular α, isIsometric_of_isIso ex43P α⟩, isLinear_of_isIso ex43P α⟩
  preStepMono {A B} φ _ := by
    refine ⟨fun {Z} g h hgh => Ex43.hom_ext ?_⟩
    have hd := congrArg (fun t : Z ⟶ B => Ex43Hom.deg t) hgh
    rw [Ex43.comp_deg, Ex43.comp_deg] at hd
    exact mul_left_cancel hd
  preStepFactor {A B} φ hφ :=
    ⟨B, φ, 𝟙 B, (Category.comp_id φ).symm, ex43_coAngular φ, hφ,
      isIsometric_of_isIso ex43P (𝟙 B), isPreStep_of_isIso ex43P (𝟙 B)⟩
  preStepFactorUniq {A B} X X' β α β' α' heq _ _ hαm hαs _ _ hα'm hα's := by
    have hαd : α.deg = 1 := hαs.1
    have hα'd : α'.deg = 1 := hα's.1
    have hd := congrArg (fun t : A ⟶ B => Ex43Hom.deg t) heq
    rw [Ex43.comp_deg, Ex43.comp_deg, hαd, hα'd] at hd
    have hd' : β.deg = β'.deg := by simpa using hd
    have hX : X = B := ex43_eq_of_isometric_preStep α hαm hαs
    have hX' : X' = B := ex43_eq_of_isometric_preStep α' hα'm hα's
    subst hX
    subst hX'
    refine ⟨Iso.refl _, Ex43.hom_ext ?_, Ex43.hom_ext ?_⟩
    · rw [hα'd, Ex43.comp_deg, hαd]
      simp
    · rw [Ex43.comp_deg]
      simp [← hd']
  preStepFactor' {A B} φ hφ :=
    ⟨A, 𝟙 A, φ, (Category.id_comp φ).symm, isIsometric_of_isIso ex43P (𝟙 A),
      isPreStep_of_isIso ex43P (𝟙 A), ex43_coAngular φ, hφ⟩
  preStepFactorUniq' {A B} X X' β α β' α' heq hβm hβs _ _ hβ'm hβ's _ _ := by
    have hβd : β.deg = 1 := hβs.1
    have hβ'd : β'.deg = 1 := hβ's.1
    have hd := congrArg (fun t : A ⟶ B => Ex43Hom.deg t) heq
    rw [Ex43.comp_deg, Ex43.comp_deg, hβd, hβ'd] at hd
    have hd' : α.deg = α'.deg := by simpa using hd
    have hX : A = X := ex43_eq_of_isometric_preStep β hβm hβs
    have hX' : A = X' := ex43_eq_of_isometric_preStep β' hβ'm hβ's
    subst hX
    subst hX'
    refine ⟨Iso.refl _, Ex43.hom_ext ?_, Ex43.hom_ext ?_⟩
    · rw [← hd']
      simp
    · rw [Ex43.comp_deg, hβ'd, hβd]
      simp
  faithfulUpToUnits {A B} φ ψ _ _ _ hφs _ hψs :=
    ⟨1, (OTimes ex43P B).one_mem, Ex43.hom_ext (by
      show φ.deg = (ψ ≫ 𝟙 B).deg
      rw [Ex43.comp_deg, Ex43.id_deg, one_mul]
      exact (show φ.deg = 1 from hφs.1).trans (show ψ.deg = 1 from hψs.1).symm)⟩
  isotropicHullExists A :=
    ⟨A, 𝟙 A, isIsometric_of_isIso ex43P (𝟙 A), isPreStep_of_isIso ex43P (𝟙 A),
      ex43_isotropicType A, fun Cc _ γ => ⟨γ, (Category.id_comp γ).symm,
        fun y hy => by rw [hy, Category.id_comp]⟩⟩
  isotropicClosed φ _ := ex43_isotropicType _

/-! ## ★10. (iii)(d) —— `𝒞^coa-pre` と `Order(Φ(A))` の 2 本の圏同値

★★`Example 4.3` では **co-angular pre-step ⟺ 次数 1**、そして
`Div ⟨1⟩ = b − a` なので、
`_A(𝒞^coa-pre)` は「`a ≤ b` なる `b`」の順序集合そのものになる。 -/

theorem ex43_coaPre_iff {a b : Ex43} (f : a ⟶ b) :
    coaPreProp ex43P f ↔ f.deg = 1 :=
  ⟨fun h => h.2.1, fun h => ⟨ex43_coAngular f, (ex43_isPreStep_iff f).mpr h⟩⟩

local instance ex43_coaPre_mult : MorphismProperty.IsMultiplicative (coaPreProp ex43P) :=
  coaPreProp_isMultiplicative ex43P ex43_core.coAngularComp

/-- ★次数 1 の射の零因子は `b − a`。 -/
theorem ex43Div_deg_one {a b : Ex43} (f : a ⟶ b) (hf : f.deg = 1) :
    ((ex43Div f : ℚ≥0) : ℚ) = b.val - a.val := by
  rw [ex43Div_val, hf]
  push_cast
  ring

/-! ### ★前置 `_A(𝒞^coa-pre) → Order(Φ(A))` -/

theorem ex43_coaPreUnder_faithful (A : Ex43) :
    (coaPreUnderFunctor ex43P A).Faithful where
  map_injective {Z W} {f g} _ := by
    refine Under.UnderMorphism.ext (InducedWideCategory.Hom.ext ?_)
    exact Ex43.hom_ext (by
      rw [(ex43_coaPre_iff _).mp f.right.property,
        (ex43_coaPre_iff _).mp g.right.property])

theorem ex43_coaPreUnder_full (A : Ex43) : (coaPreUnderFunctor ex43P A).Full := by
  constructor
  intro Z W h
  have hZ : Z.hom.hom.deg = 1 := (ex43_coaPre_iff _).mp Z.hom.property
  have hW : W.hom.hom.deg = 1 := (ex43_coaPre_iff _).mp W.hom.property
  have hle : MLe (ex43Div Z.hom.hom) (ex43Div W.hom.hom) := leOfHom h
  have hval : Z.right.obj.val ≤ W.right.obj.val := by
    obtain ⟨c, hc⟩ := hle
    have h2 : ((ex43Div Z.hom.hom : ℚ≥0) : ℚ) + (c : ℚ)
        = ((ex43Div W.hom.hom : ℚ≥0) : ℚ) := by
      rw [← NNRat.coe_add]
      exact congrArg (fun t : ℚ≥0 => (t : ℚ)) hc
    rw [ex43Div_deg_one _ hZ, ex43Div_deg_one _ hW] at h2
    have hc0 : (0 : ℚ) ≤ (c : ℚ) := c.coe_nonneg
    linarith
  obtain ⟨t, ht⟩ : ∃ t : Z.right.obj ⟶ W.right.obj, t.deg = 1 :=
    ⟨⟨1, by simpa using hval⟩, rfl⟩
  have hcomp : Z.hom.hom ≫ t = W.hom.hom :=
    Ex43.hom_ext (by rw [Ex43.comp_deg, ht, hZ, hW, mul_one])
  refine ⟨Under.homMk (show Z.right ⟶ W.right from
    ⟨t, (ex43_coaPre_iff t).mpr ht⟩) (WideSubcategory.hom_ext _ hcomp), ?_⟩
  exact Subsingleton.elim _ _

theorem ex43_coaPreUnder_essSurj (A : Ex43) : (coaPreUnderFunctor ex43P A).EssSurj := by
  refine ⟨fun x => ?_⟩
  obtain ⟨y, hy⟩ : ∃ y : ℚ≥0, toOrderCat y = x := ⟨x, rfl⟩
  obtain ⟨φ, hφ⟩ : ∃ φ : A ⟶ Ex43.mk (A.val + (y : ℚ)), φ.deg = 1 :=
    ⟨⟨1, by simpa using y.coe_nonneg⟩, rfl⟩
  refine ⟨Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp ex43P)) ⟶ ⟨_⟩ from
    ⟨φ, (ex43_coaPre_iff φ).mpr hφ⟩), ⟨eqToIso ?_⟩⟩
  rw [← hy]
  show toOrderCat (ex43Div φ) = toOrderCat y
  refine congrArg toOrderCat (NNRat.coe_injective ?_)
  rw [ex43Div_deg_one φ hφ]
  simp

theorem ex43_coaPreUnderEquiv (A : Ex43) : (coaPreUnderFunctor ex43P A).IsEquivalence :=
  ⟨ex43_coaPreUnder_faithful A, ex43_coaPreUnder_full A, ex43_coaPreUnder_essSurj A⟩

/-! ### ★後置 `(𝒞^coa-pre)_A → Order(Φ(A))^opp` -/

theorem ex43_coaPreOver_faithful (A : Ex43) :
    (coaPreOverFunctor ex43P A).Faithful where
  map_injective {Z W} {f g} _ := by
    refine Over.OverMorphism.ext (InducedWideCategory.Hom.ext ?_)
    exact Ex43.hom_ext (by
      rw [(ex43_coaPre_iff _).mp f.left.property,
        (ex43_coaPre_iff _).mp g.left.property])

/-- ★後置の関手の値は `A − B`(底は `𝟙` なので `Φ.map` は恒等)。 -/
theorem ex43_coaPreOver_obj {A : Ex43} (Z : Over (⟨A⟩ : WideSubcategory (coaPreProp ex43P))) :
    ((coaPreOverFunctor ex43P A).obj Z).unop = toOrderCat (ex43Div Z.hom.hom) := rfl

theorem ex43_coaPreOver_full (A : Ex43) : (coaPreOverFunctor ex43P A).Full := by
  constructor
  intro Z W h
  have hZ : Z.hom.hom.deg = 1 := (ex43_coaPre_iff _).mp Z.hom.property
  have hW : W.hom.hom.deg = 1 := (ex43_coaPre_iff _).mp W.hom.property
  have hle : MLe (ex43Div W.hom.hom) (ex43Div Z.hom.hom) := leOfHom h.unop
  have hval : Z.left.obj.val ≤ W.left.obj.val := by
    obtain ⟨c, hc⟩ := hle
    have h2 : ((ex43Div W.hom.hom : ℚ≥0) : ℚ) + (c : ℚ)
        = ((ex43Div Z.hom.hom : ℚ≥0) : ℚ) := by
      rw [← NNRat.coe_add]
      exact congrArg (fun t : ℚ≥0 => (t : ℚ)) hc
    rw [ex43Div_deg_one _ hZ, ex43Div_deg_one _ hW] at h2
    have hc0 : (0 : ℚ) ≤ (c : ℚ) := c.coe_nonneg
    linarith
  obtain ⟨t, ht⟩ : ∃ t : Z.left.obj ⟶ W.left.obj, t.deg = 1 :=
    ⟨⟨1, by simpa using hval⟩, rfl⟩
  have hcomp : t ≫ W.hom.hom = Z.hom.hom :=
    Ex43.hom_ext (by rw [Ex43.comp_deg, ht, hZ, hW, one_mul])
  refine ⟨Over.homMk (show Z.left ⟶ W.left from
    ⟨t, (ex43_coaPre_iff t).mpr ht⟩) (WideSubcategory.hom_ext _ hcomp), ?_⟩
  exact Subsingleton.elim _ _

theorem ex43_coaPreOver_essSurj (A : Ex43) : (coaPreOverFunctor ex43P A).EssSurj := by
  refine ⟨fun x => ?_⟩
  obtain ⟨y, hy⟩ : ∃ y : ℚ≥0, Opposite.op (toOrderCat y) = x := ⟨x.unop, rfl⟩
  obtain ⟨φ, hφ⟩ : ∃ φ : Ex43.mk (A.val - (y : ℚ)) ⟶ A, φ.deg = 1 :=
    ⟨⟨1, by simpa using y.coe_nonneg⟩, rfl⟩
  refine ⟨Over.mk (show (⟨_⟩ : WideSubcategory (coaPreProp ex43P)) ⟶ ⟨A⟩ from
    ⟨φ, (ex43_coaPre_iff φ).mpr hφ⟩), ⟨eqToIso ?_⟩⟩
  rw [← hy]
  refine Opposite.unop_injective ?_
  show toOrderCat (ex43Div φ) = toOrderCat y
  refine congrArg toOrderCat (NNRat.coe_injective ?_)
  rw [ex43Div_deg_one φ hφ]
  simp

theorem ex43_coaPreOverEquiv (A : Ex43) : (coaPreOverFunctor ex43P A).IsEquivalence :=
  ⟨ex43_coaPreOver_faithful A, ex43_coaPreOver_full A, ex43_coaPreOver_essSurj A⟩

/-- ★★★★**`Example 4.3` の `𝒞` は Frobenioid**。 -/
theorem ex43_frobenioid : Frobenioid ex43P where
  core := ex43_core
  coaPreUnderEquiv := ex43_coaPreUnderEquiv
  coaPreOverEquiv := ex43_coaPreOverEquiv

/-! ## ★11. group-like でない・底圏と `Φ` の型

原文 (FrdI p.82):
> no object of C is group-like. Thus, one
-/

/-- ★`ℚ≥0` は **group-like でない** —— sharp なので、`1` が可逆なら `1 = 0`。 -/
theorem not_isGroupLike_nnrat : ¬ IsGroupLike ℚ≥0 := by
  intro h
  have h1 : IsAddUnit (1 : ℚ≥0) := (isGroupLike_iff ℚ≥0).mp h 1
  exact one_ne_zero (isSharp_nnrat 1 h1)

/-- ★★**主張 8** —— `𝒞` の対象はどれも **group-like でない**。 -/
theorem ex43_not_groupLikeType : ¬ IsOfGroupLikeType ex43P := fun h =>
  not_isGroupLike_nnrat (h (Ex43.mk 0))

/-- ★底圏 `𝒟` の射はすべて同型(`Discrete` の射は等式)。 -/
theorem d43_isIso {A B : Discrete PUnit} (f : A ⟶ B) : IsIso f := by
  obtain rfl : A = B := Subsingleton.elim _ _
  exact ⟨f, rfl, rfl⟩

/-- ★★**主張 10** —— `𝒟` は **FSM-type**。 -/
theorem d43_fsmType : IsOfFSMType (Discrete PUnit) := fun _ _ β _ => d43_isIso β

/-- ★★**主張 10** —— したがって **FSMFF-type**。

原文 (FrdI p.82):
> FSMFF-type, and
-/
theorem d43_fsmffType : IsOfFSMFFType (Discrete PUnit) :=
  isOfFSMFFType_of_isOfFSMType d43_fsmType

/-- ★★**主張 13** —— `𝒟` は **slim**(射が 1 本しかないのでスライス圏の自己自然変換は恒等)。

原文 (FrdI p.82):
> a slim base category D. Now one verifies immediately that if
-/
theorem d43_slim : IsSlimCat (Discrete PUnit) := by
  intro A η
  refine Iso.ext (NatTrans.ext (funext fun X => ?_))
  exact Subsingleton.elim _ _

/-- ★★**主張 11** —— `Φ` は **non-dilating**(定数関手なので誘導射は恒等)。 -/
theorem phi43_nonDilating : MonoidOn.IsNonDilatingOn Phi43 := by
  intro A e
  have h : Phi43.map e = AddMonoidHom.id ℚ≥0 := by ext x; rfl
  rw [h]
  exact isNonDilating_id (M := ℚ≥0)

/-! ## ★12. Frobenius-normalized 型と standard 型 -/

/-- ★`End a` の冪の次数は次数の冪(`End` の積は `x * y = y ≫ x`)。 -/
theorem ex43_end_pow_deg {a : Ex43} (α : End a) (n : ℕ) :
    (α ^ n : End a).deg = α.deg ^ n := by
  induction n with
  | zero => rfl
  | succ k ih =>
      have h1 : ((α ^ (k + 1) : End a)).deg = ((α ≫ (α ^ k : End a) : a ⟶ a)).deg :=
        congrArg Ex43Hom.deg (pow_succ α k)
      rw [h1, Ex43.comp_deg, ih, pow_succ]

/-- ★★`𝒞` は **Frobenius-normalized 型** —— `𝒪^▷` が自明だから。 -/
theorem ex43_frobNormalizedType : IsOfFrobeniusNormalizedType ex43P := by
  intro A φ _ α hα
  have hd : α.deg = 1 := hα.2
  refine Ex43.hom_ext ?_
  rw [Ex43.comp_deg, Ex43.comp_deg, ex43_end_pow_deg, hd, one_pow, one_mul, mul_one]

/-- ★★★**主張 12** —— `𝒞` は **standard 型**。

★(a) isotropic 型なので quasi-isotropic 型、かつ恒等射が Frobenius 型なので
Frobenius-isotropic 型。★(b) group-like 型でないので前件が偽。
★(c)(d)(e) は上で取った。

原文 (FrdI p.82):
> it follows that C is also of standard type, over
-/
theorem ex43_standardType : IsOfStandardType (Discrete PUnit) Ex43 ex43P ex43_core where
  quasiIsotropic :=
    isOfQuasiIsotropicType_of_isOfIsotropicType ex43P ex43_core ex43_isotropicType
  frobIsotropic := fun A =>
    ⟨A, 𝟙 A, isFrobeniusType_of_isIso ex43P (𝟙 A), ex43_isotropicType A⟩
  groupLikeCompact := fun hgl => absurd hgl ex43_not_groupLikeType
  frobNormalized := ex43_frobNormalizedType
  baseFSMFF := d43_fsmffType
  phiNonDilating := phi43_nonDilating

/-! ## ★13. 自己同値 `Ψ_λ`

原文 (FrdI p.82):
> determines a self-equivalence of categories

★**`a ↦ a`(`a ≥ 0`)・`−a ↦ −λ·a`** を 1 本の式にまとめると
`ψ_λ(q) = if 0 ≤ q then q else λ · q` である。 -/

/-- ★`Ψ_λ` の対象上の作用。 -/
def psiVal (lam : ℚ) (q : ℚ) : ℚ := if 0 ≤ q then q else lam * q

@[simp] theorem psiVal_of_nonneg {lam q : ℚ} (h : 0 ≤ q) : psiVal lam q = q := if_pos h

@[simp] theorem psiVal_of_neg {lam q : ℚ} (h : q < 0) : psiVal lam q = lam * q :=
  if_neg (not_le.mpr h)

/-- ★`ψ_λ` は符号を保つ。 -/
theorem psiVal_neg_of_neg {lam q : ℚ} (hl : 0 < lam) (h : q < 0) : psiVal lam q < 0 := by
  rw [psiVal_of_neg h]
  exact mul_neg_of_pos_of_neg hl h

/-- ★★`ψ_λ` は「`n · a ≤ b`」を保つ。 -/
theorem psiVal_cond {lam : ℚ} (hl : 0 < lam) {n a b : ℚ} (hn : 0 < n)
    (h : n * a ≤ b) : n * psiVal lam a ≤ psiVal lam b := by
  rcases le_or_gt 0 a with ha | ha
  · rcases le_or_gt 0 b with hb | hb
    · rwa [psiVal_of_nonneg ha, psiVal_of_nonneg hb]
    · exact absurd h (not_le.mpr (lt_of_lt_of_le hb (mul_nonneg hn.le ha)))
  · rcases le_or_gt 0 b with hb | hb
    · rw [psiVal_of_nonneg hb, psiVal_of_neg ha]
      exact le_trans (mul_neg_of_pos_of_neg hn (mul_neg_of_pos_of_neg hl ha)).le hb
    · rw [psiVal_of_neg ha, psiVal_of_neg hb]
      calc n * (lam * a) = lam * (n * a) := by ring
        _ ≤ lam * b := mul_le_mul_of_nonneg_left h hl.le

/-- ★★`ψ_λ` は「`n · a ≤ b`」を**反射する**。 -/
theorem psiVal_cond_rev {lam : ℚ} (hl : 0 < lam) {n a b : ℚ} (hn : 0 < n)
    (h : n * psiVal lam a ≤ psiVal lam b) : n * a ≤ b := by
  rcases le_or_gt 0 a with ha | ha
  · rcases le_or_gt 0 b with hb | hb
    · rwa [psiVal_of_nonneg ha, psiVal_of_nonneg hb] at h
    · rw [psiVal_of_nonneg ha, psiVal_of_neg hb] at h
      exact absurd h (not_le.mpr (lt_of_lt_of_le (mul_neg_of_pos_of_neg hl hb)
        (mul_nonneg hn.le ha)))
  · rcases le_or_gt 0 b with hb | hb
    · exact le_trans (mul_neg_of_pos_of_neg hn ha).le hb
    · rw [psiVal_of_neg ha, psiVal_of_neg hb] at h
      have h2 : lam * (n * a) ≤ lam * b := by linarith [h]
      exact le_of_mul_le_mul_left h2 hl

/-- ★★★**自己関手 `Ψ_λ`** —— 対象は `ψ_λ`、射は次数をそのまま写す。 -/
def Psi43 (lam : ℚ) (hl : 0 < lam) : Ex43 ⥤ Ex43 where
  obj a := Ex43.mk (psiVal lam a.val)
  map {a b} f := ⟨f.deg, psiVal_cond hl (by exact_mod_cast f.deg.pos) f.cond⟩
  map_id a := Ex43.hom_ext rfl
  map_comp f g := Ex43.hom_ext rfl

/-- ★★**`Ψ_λ` は Frobenius 次数を保つ**。

原文 (FrdI p.82):
> that preserves Frobenius degrees [cf. Theorem 3.4, (iii)]. On the other hand, it
-/
@[simp] theorem Psi43_map_deg {lam : ℚ} (hl : 0 < lam) {a b : Ex43} (f : a ⟶ b) :
    ((Psi43 lam hl).map f).deg = f.deg := rfl

theorem Psi43_degFr {lam : ℚ} (hl : 0 < lam) {a b : Ex43} (f : a ⟶ b) :
    ex43P.degFr ((Psi43 lam hl).map f) = ex43P.degFr f := rfl

theorem Psi43_faithful {lam : ℚ} (hl : 0 < lam) : (Psi43 lam hl).Faithful where
  map_injective {a b} {f g} h := Ex43.hom_ext (by
    have h2 : ((Psi43 lam hl).map f).deg = ((Psi43 lam hl).map g).deg :=
      congrArg Ex43Hom.deg h
    exact h2)

theorem Psi43_full {lam : ℚ} (hl : 0 < lam) : (Psi43 lam hl).Full := by
  constructor
  intro a b g
  refine ⟨⟨g.deg, psiVal_cond_rev hl (by exact_mod_cast g.deg.pos) g.cond⟩, ?_⟩
  exact Ex43.hom_ext rfl

theorem Psi43_essSurj {lam : ℚ} (hl : 0 < lam) : (Psi43 lam hl).EssSurj := by
  refine ⟨fun q => ?_⟩
  refine ⟨Ex43.mk (if 0 ≤ q.val then q.val else lam⁻¹ * q.val), ⟨eqToIso (Ex43.ext ?_)⟩⟩
  show psiVal lam (if 0 ≤ q.val then q.val else lam⁻¹ * q.val) = q.val
  rcases le_or_gt 0 q.val with h | h
  · rw [if_pos h, psiVal_of_nonneg h]
  · rw [if_neg (not_le.mpr h)]
    have hneg : lam⁻¹ * q.val < 0 := mul_neg_of_pos_of_neg (inv_pos.mpr hl) h
    rw [psiVal_of_neg hneg, ← mul_assoc, mul_inv_cancel₀ hl.ne', one_mul]

/-- ★★★**`Ψ_λ` は自己同値**。 -/
theorem Psi43_isEquivalence {lam : ℚ} (hl : 0 < lam) : (Psi43 lam hl).IsEquivalence :=
  ⟨Psi43_faithful hl, Psi43_full hl, Psi43_essSurj hl⟩

/-! ## ★14. `Ψ_λ` は `λ ≠ 1` なら恒等関手と同型でない

★★これが原文の「**右辺の同型は恒等・左辺の同型は `λ` 倍**」の内容である。
`Ψ_λ` は Frobenius 次数を保つ自己同値でありながら、`λ ≠ 1` では恒等と同型でない。

原文 (FrdI p.82):
> of Theorem 4.2, (iii), is the identity on
-/

/-- ★同型射があれば対象は等しい(同型は等長 pre-step)。 -/
theorem ex43_eq_of_isIso {a b : Ex43} (f : a ⟶ b) [IsIso f] : a = b :=
  ex43_eq_of_isometric_preStep f (isIsometric_of_isIso ex43P f) (isPreStep_of_isIso ex43P f)

/-- ★★★`λ ≠ 1` なら `Ψ_λ ≇ 𝟭` —— `−1` の像が `−λ` だから。 -/
theorem Psi43_not_iso_id {lam : ℚ} (hl : 0 < lam) (hne : lam ≠ 1) :
    IsEmpty (Psi43 lam hl ≅ 𝟭 Ex43) := by
  constructor
  intro η
  have hobj : (Psi43 lam hl).obj (Ex43.mk (-1)) = Ex43.mk (-1) := by
    have hi := (η.app (Ex43.mk (-1))).isIso_hom
    exact ex43_eq_of_isIso (η.app (Ex43.mk (-1))).hom
  have hval : psiVal lam (-1 : ℚ) = (-1 : ℚ) := congrArg Ex43.val hobj
  rw [psiVal_of_neg (by norm_num : (-1 : ℚ) < 0)] at hval
  exact hne (by linarith)

/-! ## ★★★★`Example 4.3` —— まとめ -/

def ex43_frobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 81, item := "Example 4.3",
    sectionId := "frdi-example-4-3" }

end ABC3.Found.FrdI
