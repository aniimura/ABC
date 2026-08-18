import ABC3.Found.FrdI.Prop44Pre

/-!
# [FrdI] Proposition 4.4 の図式の中段 —— `𝒞^birat → 𝔽_{Φ^gp}` の零因子

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.83。

原文 (FrdI p.83):
> where the functors between elementary Frobenioids are those induced by the natural

★`Prop44.lean` は図式の**下段** `𝒞^birat → 𝔽_{0_𝒟}` までを作っている。
★★**中段 `𝒞^birat → 𝔽_{Φ^gp}` が無い**ので、`Proposition 4.4, (iii)` の
`Φ^birat`(その像として定まる)も書けない —— ★**ここが §4 の律速である**。

## ★★中段の零因子の式(紙の上で確定、2026-08-17)

★代表元 `(a : A′ → A, φ : A′ → B)` は形式的な合成 `a⁻¹ ≫ φ` を表す。
`𝔽_{Φ^gp}` の合成則(`Div_comp`)から

```
Div(a⁻¹)      = − Φ^gp(Base a)⁻¹ (Div a)
Div(a⁻¹ ≫ φ) = Φ^gp(Base a)⁻¹ (Div φ) − degFr(φ) · Φ^gp(Base a)⁻¹ (Div a)
```

★★**`degFr φ` の係数を落とさないこと** —— `Div_comp` の第 2 項が
`degFr φ • Div ψ` だからである(初稿でここを落としかけた)。

## ★★代表元の取り替えで不変であること

★遷移 `c`(co-angular pre-step、`c ≫ a_Z = a_W`)で `φ ↦ c ≫ φ` としたとき、
`Div c` の項が **2 か所に現れて打ち消す** ——
第 1 項の `degFr φ • Div c` と、第 2 項の `degFr φ · (Div a_W の Div c 成分)` である。
★★これが「不変量である」ことの中身で、**`a` が pre-step(次数 1)**が効く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★0. `toGp` は加法的 -/

variable {M : Type w} [AddCommMonoid M]

@[simp] theorem toGp_zero : toGp M 0 = 0 := by
  show (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) 0 = 0
  exact map_zero _

theorem toGp_add (x y : M) : toGp M (x + y) = toGp M x + toGp M y := by
  show (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) (x + y)
    = (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) x
      + (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) y
  exact map_add _ x y

theorem toGp_nsmul (n : ℕ) (x : M) : toGp M (n • x) = n • toGp M x := by
  induction n with
  | zero => simp
  | succ k ih => rw [succ_nsmul, toGp_add, ih, succ_nsmul]

/-! ## ★1. 代表元が定める `Φ^gp(A_𝒟)` の元 -/

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

variable {P} in
/-- ★★**代表元 `(a, φ)` が定める `Φ^gp(A_𝒟)` の元** —— 形式的な合成 `a⁻¹ ≫ φ` の零因子。

★`sliceBaseOf` と同じく **`IsIso` を仮引数に割る**(`Prop44.lean` の測定)。 -/
noncomputable def sliceDivGpOf {A B A' : C} (a : A' ⟶ A) (_ha : IsIso (P.Base a))
    (φ : A' ⟶ B) : Gp (Φ.val (P.toElem.obj A).base) :=
  haveI := _ha
  gpMap _ (Φ.map (inv (P.Base a)))
    (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div a))

variable {P} in
theorem sliceDivGpOf_eq {A B A' : C} (a : A' ⟶ A) (ha : IsIso (P.Base a)) (φ : A' ⟶ B) :
    sliceDivGpOf a ha φ = haveI := ha
      gpMap _ (Φ.map (inv (P.Base a)))
        (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div a)) := rfl

variable {P} in
/-- ★`IsIso` は `Prop` なので、射が等しければ値も等しい。 -/
theorem sliceDivGpOf_congr {A B A' : C} {a a' : A' ⟶ A} (h : a = a') (ha : IsIso (P.Base a))
    (ha' : IsIso (P.Base a')) (φ : A' ⟶ B) :
    sliceDivGpOf a ha φ = sliceDivGpOf a' ha' φ := by
  subst h; rfl

/-! ## ★2. ★★★代表元の取り替えで不変 -/

variable (Φ) in
/-- ★`gpMap ∘ gpMap` を `Φ.map` の合成にまとめる。 -/
theorem gpMap_phi_comp {A B E : D} (β : A ⟶ B) (α : B ⟶ E) (x : Gp (Φ.val E)) :
    gpMap _ (Φ.map β) (gpMap _ (Φ.map α) x) = gpMap _ (Φ.map (β ≫ α)) x := by
  have h : (Φ.map β).comp (Φ.map α) = Φ.map (β ≫ α) := by
    ext y
    exact (Φ.map_comp α β y).symm
  rw [← h, gpMap_comp]
  rfl

variable {P} in
/-- ★★★**不変量である** —— 遷移 `c` で代表元を取り替えても値が変わらない。

★★`Div c` の項が**第 1 項と第 2 項の両方に現れて打ち消す**のが要点。 -/
theorem sliceDivGpOf_precomp {A B A' A'' : C} (a : A' ⟶ A) (ha : IsIso (P.Base a))
    (has : P.degFr a = 1) (c : A'' ⟶ A') (_hc : IsIso (P.Base c))
    (hca : IsIso (P.Base (c ≫ a))) (φ : A' ⟶ B) (hcs : P.degFr c = 1) :
    sliceDivGpOf (c ≫ a) hca (c ≫ φ) = sliceDivGpOf a ha φ := by
  haveI := ha
  haveI := hca
  -- ★分子を `Φ.map (Base c)` の像としてまとめる
  have hnum : toGp _ (P.Div (c ≫ φ))
        - ((P.degFr (c ≫ φ) : ℕ+) : ℕ) • toGp _ (P.Div (c ≫ a))
      = gpMap _ (Φ.map (P.Base c))
        (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div a)) := by
    rw [P.Div_comp, P.Div_comp, P.degFr_comp, hcs, mul_one, has, PNat.one_coe, one_smul,
      toGp_add, toGp_add, toGp_nsmul, map_sub, gpMap_toGp, map_nsmul, gpMap_toGp]
    simp only [smul_add]
    abel
  have hkey : inv (P.Base (c ≫ a)) ≫ P.Base c = inv (P.Base a) := by
    rw [IsIso.inv_comp_eq, P.Base_comp, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
  rw [sliceDivGpOf_eq, sliceDivGpOf_eq, hnum, gpMap_phi_comp]
  exact congrArg (fun t => gpMap _ (Φ.map t)
    (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div a))) hkey

/-! ## ★3. `Hom^birat` の上の関数として実現する

★★**測定(2026-08-17)**: 余錐(`Cocone`)では作れない ——
`homFunctorBirat` の値は `Type (max v u2 v2)` に住むのに、
`Gp (Φ.val −)` は **`Type w`** に住み、`w` は他と独立だからである。
★`biratDeg`(`ℕ+`)や `biratBase`(`Type v`)は `ULift` で持ち上がるが、
★★**`Div` だけは持ち上がらない**。

★★★**逃げ道: 代表元を選んで定義し、「代表元によらない」を別に示す。**
`HomBirat.eq_iff` が共通の上界での一致を与えるので、
不変性(`sliceDivGpOf_precomp`)を 2 回当てれば済む。 -/

variable (G : Frobenioid P)

variable {P G} in
/-- ★★**代表元によらない** —— 共通の上界へ両方を送って不変性を当てる。 -/
theorem sliceDivGpOf_eq_of_mk_eq {A B : C} {Z W : IdxBirat P G A}
    {φ : Z.unop.left.obj ⟶ B} {ψ : W.unop.left.obj ⟶ B}
    (h : HomBirat.mk Z φ = HomBirat.mk W ψ) :
    sliceDivGpOf (P := P) Z.unop.hom.hom Z.unop.hom.property.2.2 φ
      = sliceDivGpOf (P := P) W.unop.hom.hom W.unop.hom.property.2.2 ψ := by
  obtain ⟨V, u, v, huv⟩ := HomBirat.eq_iff.mp h
  have hZ : u.unop.left.hom ≫ Z.unop.hom.hom = V.unop.hom.hom :=
    congrArg (fun t : V.unop.left ⟶ coaPreObj P G A => t.hom) (Over.w u.unop)
  have hW : v.unop.left.hom ≫ W.unop.hom.hom = V.unop.hom.hom :=
    congrArg (fun t : V.unop.left ⟶ coaPreObj P G A => t.hom) (Over.w v.unop)
  have hcaZ : IsIso (P.Base (u.unop.left.hom ≫ Z.unop.hom.hom)) := by
    rw [hZ]; exact V.unop.hom.property.2.2
  have hcaW : IsIso (P.Base (v.unop.left.hom ≫ W.unop.hom.hom)) := by
    rw [hW]; exact V.unop.hom.property.2.2
  have h1 : sliceDivGpOf (P := P) V.unop.hom.hom V.unop.hom.property.2.2
        (u.unop.left.hom ≫ φ)
      = sliceDivGpOf (P := P) Z.unop.hom.hom Z.unop.hom.property.2.2 φ :=
    (sliceDivGpOf_congr hZ hcaZ V.unop.hom.property.2.2 (u.unop.left.hom ≫ φ)).symm.trans
      (sliceDivGpOf_precomp Z.unop.hom.hom Z.unop.hom.property.2.2
        Z.unop.hom.property.2.1 u.unop.left.hom u.unop.left.property.2.2 hcaZ φ
        u.unop.left.property.2.1)
  have h2 : sliceDivGpOf (P := P) V.unop.hom.hom V.unop.hom.property.2.2
        (v.unop.left.hom ≫ ψ)
      = sliceDivGpOf (P := P) W.unop.hom.hom W.unop.hom.property.2.2 ψ :=
    (sliceDivGpOf_congr hW hcaW V.unop.hom.property.2.2 (v.unop.left.hom ≫ ψ)).symm.trans
      (sliceDivGpOf_precomp W.unop.hom.hom W.unop.hom.property.2.2
        W.unop.hom.property.2.1 v.unop.left.hom v.unop.left.property.2.2 hcaW ψ
        v.unop.left.property.2.1)
  rw [← h1, ← h2, huv]

variable {P G} in
/-- ★★★**`Hom^birat` の零因子**(`Φ^gp` の中で)—— 原文の図式の中段。 -/
noncomputable def biratDivGp {A B : C} (f : HomBirat P G A B) :
    Gp (Φ.val (P.toElem.obj A).base) :=
  sliceDivGpOf (P := P) (HomBirat.exists_rep f).choose.unop.hom.hom
    (HomBirat.exists_rep f).choose.unop.hom.property.2.2
    (HomBirat.exists_rep f).choose_spec.choose

variable {P G} in
@[simp] theorem biratDivGp_mk {A B : C} (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) :
    biratDivGp (HomBirat.mk Z φ)
      = sliceDivGpOf (P := P) Z.unop.hom.hom Z.unop.hom.property.2.2 φ :=
  sliceDivGpOf_eq_of_mk_eq (HomBirat.exists_rep (HomBirat.mk Z φ)).choose_spec.choose_spec

variable {P G} in
/-- ★**`𝒞` から来た射の零因子はもとのまま**(`𝔽_Φ → 𝔽_{Φ^gp}` の像)。 -/
@[simp] theorem biratDivGp_toHomBirat {A B : C} (φ : A ⟶ B) :
    biratDivGp (toHomBirat (P := P) (G := G) φ) = toGp _ (P.Div φ) := by
  rw [toHomBirat, biratDivGp_mk]
  show sliceDivGpOf (P := P) (𝟙 A) _ φ = _
  haveI : IsIso (P.Base (𝟙 A)) := by rw [P.Base_id]; infer_instance
  rw [sliceDivGpOf_eq]
  show gpMap _ (Φ.map (inv (P.Base (𝟙 A))))
      (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div (𝟙 A))) = _
  have hd : P.Div (𝟙 A) = 0 := P.Div_id A
  have hb : inv (P.Base (𝟙 A)) = 𝟙 _ := by
    refine (IsIso.inv_eq_of_hom_inv_id ?_)
    rw [P.Base_id, Category.id_comp]
  rw [hd, hb]
  have hmap : Φ.map (𝟙 ((P.toElem.obj A).base)) = AddMonoidHom.id _ := by
    ext y; exact Φ.map_id _ y
  rw [hmap]
  simp [gpMap_id]

/-! ## ★4. ★★★★合成則

★★**紙の上の計算(2026-08-17)**。`f = (a, φ)`、`g = (b, ψ)`、引き戻しを `(γ, α)`
(`γ ≫ φ = α ≫ b`)とし、`m = degFr φ`、`n = degFr ψ` と置く。

★**両辺とも `Φ^gp(Base a)⁻¹` の像**なので、`A′` の上の等式に落ちる。
★さらに `Φ^gp(Base γ)⁻¹` を外に出すと、**`Dd` の上の等式**

```
Div(α ≫ ψ) − (n·m)·Div(γ ≫ a) = Φ^gp(Base α)(…ψ…) + n · Φ^gp(Base γ)(…φ…)
```

になり、★★**残るのは `γ ≫ φ = α ≫ b` の零因子成分**

```
Φ(Base γ)(Div φ) + m·Div γ = Φ(Base α)(Div b) + Div α
```

**だけ**である。★これが「合成が well-defined」の中身そのものである。 -/

variable {P G} in
/-- ★引き戻しの `α` の次数は `φ` の次数に等しい —— `γ`・`b` が pre-step だから。 -/
theorem biratPullAlpha_degFr' (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    P.degFr (biratPullAlpha F Z φ W) = P.degFr φ := by
  have h : P.degFr φ * P.degFr (biratPullGamma F Z φ W)
      = P.degFr W.unop.hom.hom * P.degFr (biratPullAlpha F Z φ W) :=
    (P.degFr_comp _ _).symm.trans
      ((congrArg P.degFr (biratPull_sq F Z φ W)).trans (P.degFr_comp _ _))
  rw [show P.degFr (biratPullGamma F Z φ W) = 1 from (biratPullGamma_preStep F Z φ W).1,
    show P.degFr W.unop.hom.hom = 1 from W.unop.hom.property.2.1, mul_one, one_mul] at h
  exact h.symm

variable {P G} in
/-- ★★**引き戻しの四角形の零因子成分** —— 合成則の核心。 -/
theorem biratPull_div (F : FrobenioidCore P) {A B : C}
    (Z : IdxBirat P G A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P G B) :
    Φ.map (P.Base (biratPullGamma F Z φ W)) (P.Div φ)
        + ((P.degFr φ : ℕ+) : ℕ) • P.Div (biratPullGamma F Z φ W)
      = Φ.map (P.Base (biratPullAlpha F Z φ W)) (P.Div W.unop.hom.hom)
        + P.Div (biratPullAlpha F Z φ W) := by
  have h : Φ.map (P.Base (biratPullGamma F Z φ W)) (P.Div φ)
        + ((P.degFr φ : ℕ+) : ℕ) • P.Div (biratPullGamma F Z φ W)
      = Φ.map (P.Base (biratPullAlpha F Z φ W)) (P.Div W.unop.hom.hom)
        + ((P.degFr W.unop.hom.hom : ℕ+) : ℕ) • P.Div (biratPullAlpha F Z φ W) :=
    (P.Div_comp _ _).symm.trans
      ((congrArg P.Div (biratPull_sq F Z φ W)).trans (P.Div_comp _ _))
  rw [show P.degFr W.unop.hom.hom = 1 from W.unop.hom.property.2.1,
    PNat.one_coe, one_smul] at h
  exact h

variable {P} in
/-- ★★★★**代表元レベルの合成則** —— `𝒞^birat → 𝔽_{Φ^gp}` の関手性の本体。

★★**両辺とも `Φ^gp(Base (γ≫a))⁻¹` の像**に書けるので、
`Dd` の上の等式に落ちる。★そこで残るのは
**`γ ≫ φ = α ≫ b` の零因子成分**だけである。 -/
theorem sliceDivGpOf_comp {A B E A' B' Dd : C}
    (a : A' ⟶ A) (ha : IsIso (P.Base a)) (has : P.degFr a = 1)
    (b : B' ⟶ B) (hb : IsIso (P.Base b)) (hbs : P.degFr b = 1)
    (φ : A' ⟶ B) (ψ : B' ⟶ E)
    (γ : Dd ⟶ A') (hγ : IsIso (P.Base γ)) (hγs : P.degFr γ = 1)
    (α : Dd ⟶ B') (hsq : γ ≫ φ = α ≫ b)
    (hca : IsIso (P.Base (γ ≫ a))) :
    sliceDivGpOf (γ ≫ a) hca (α ≫ ψ)
      = gpMap _ (Φ.map (sliceBaseOf (P := P) a ha φ)) (sliceDivGpOf b hb ψ)
        + ((P.degFr ψ : ℕ+) : ℕ) • sliceDivGpOf a ha φ := by
  haveI := ha; haveI := hb; haveI := hγ; haveI := hca
  -- ★次数
  have hdα : P.degFr α = P.degFr φ := by
    have h : P.degFr φ * P.degFr γ = P.degFr b * P.degFr α :=
      (P.degFr_comp _ _).symm.trans ((congrArg P.degFr hsq).trans (P.degFr_comp _ _))
    rw [hγs, hbs, mul_one, one_mul] at h
    exact h.symm
  -- ★底の 2 本
  have hbsq : P.Base γ ≫ P.Base φ = P.Base α ≫ P.Base b := by
    rw [← P.Base_comp, ← P.Base_comp, hsq]
  have hb1 : inv (P.Base γ) ≫ P.Base α = P.Base φ ≫ inv (P.Base b) := by
    rw [IsIso.inv_comp_eq, ← Category.assoc, hbsq, Category.assoc, IsIso.hom_inv_id,
      Category.comp_id]
  have hinvca : inv (P.Base (γ ≫ a)) = inv (P.Base a) ≫ inv (P.Base γ) := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [P.Base_comp, Category.assoc, ← Category.assoc (P.Base a), IsIso.hom_inv_id,
      Category.id_comp, IsIso.hom_inv_id]
  -- ★零因子の等式(`hsq` の成分)
  have hdiv : Φ.map (P.Base γ) (P.Div φ) + ((P.degFr φ : ℕ+) : ℕ) • P.Div γ
      = Φ.map (P.Base α) (P.Div b) + P.Div α := by
    have h : Φ.map (P.Base γ) (P.Div φ) + ((P.degFr φ : ℕ+) : ℕ) • P.Div γ
        = Φ.map (P.Base α) (P.Div b) + ((P.degFr b : ℕ+) : ℕ) • P.Div α :=
      (P.Div_comp _ _).symm.trans ((congrArg P.Div hsq).trans (P.Div_comp _ _))
    rwa [hbs, PNat.one_coe, one_smul] at h
  have hB : toGp _ (P.Div α)
      = toGp _ (Φ.map (P.Base γ) (P.Div φ)) + ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div γ)
        - toGp _ (Φ.map (P.Base α) (P.Div b)) := by
    have hgp : toGp (Φ.val (P.toElem.obj Dd).base) (Φ.map (P.Base γ) (P.Div φ))
          + ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div γ)
        = toGp _ (Φ.map (P.Base α) (P.Div b)) + toGp _ (P.Div α) := by
      rw [← toGp_nsmul, ← toGp_add, ← toGp_add, hdiv]
    exact eq_sub_of_add_eq' hgp.symm
  -- ★★★`Dd` の上の等式
  have hkey : toGp _ (P.Div (α ≫ ψ))
        - ((P.degFr (α ≫ ψ) : ℕ+) : ℕ) • toGp _ (P.Div (γ ≫ a))
      = gpMap _ (Φ.map (P.Base α))
          (toGp _ (P.Div ψ) - ((P.degFr ψ : ℕ+) : ℕ) • toGp _ (P.Div b))
        + ((P.degFr ψ : ℕ+) : ℕ) • gpMap _ (Φ.map (P.Base γ))
          (toGp _ (P.Div φ) - ((P.degFr φ : ℕ+) : ℕ) • toGp _ (P.Div a)) := by
    rw [P.Div_comp, P.Div_comp, P.degFr_comp, hdα, has, PNat.one_coe, one_smul,
      toGp_add, toGp_add, toGp_nsmul, map_sub, map_sub, gpMap_toGp, gpMap_toGp,
      map_nsmul, map_nsmul, gpMap_toGp, gpMap_toGp, hB, PNat.mul_coe]
    simp only [smul_sub, smul_add, smul_smul]
    abel
  -- ★★両辺を `Φ^gp(Base (γ≫a))⁻¹` の像に揃える
  have e1 : (inv (P.Base a) ≫ inv (P.Base γ)) ≫ P.Base α
      = (inv (P.Base a) ≫ P.Base φ) ≫ inv (P.Base b) := by
    rw [Category.assoc, hb1, ← Category.assoc]
  have e2 : (inv (P.Base a) ≫ inv (P.Base γ)) ≫ P.Base γ = inv (P.Base a) := by
    rw [Category.assoc, IsIso.inv_hom_id, Category.comp_id]
  rw [sliceDivGpOf_eq, sliceDivGpOf_eq, sliceDivGpOf_eq, sliceBaseOf_eq, hkey, hinvca,
    map_add, map_nsmul, gpMap_phi_comp, gpMap_phi_comp, gpMap_phi_comp, e1, e2]

variable {P G} in
/-- ★★★★**`Hom^birat` の零因子は合成を保つ** —— 中段の関手性の本体。

原文 (FrdI p.83):
> where the functors between elementary Frobenioids are those induced by the natural

★`𝔽_{Φ^gp}` の合成則そのものの形
`Div(f ≫ g) = Φ^gp(Base f)(Div g) + degFr(g) · Div(f)` になっている。 -/
theorem biratDivGp_comp {A B E : C} (f : HomBirat P G A B) (g : HomBirat P G B E) :
    biratDivGp (compBirat P G G.core f g)
      = gpMap _ (Φ.map (biratBase f)) (biratDivGp g)
        + ((biratDeg g : ℕ+) : ℕ) • biratDivGp f := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  rw [compBirat_mk, biratDivGp_mk, biratDivGp_mk, biratDivGp_mk, biratBase_mk, biratDeg_mk]
  exact sliceDivGpOf_comp Z.unop.hom.hom Z.unop.hom.property.2.2 Z.unop.hom.property.2.1
    W.unop.hom.hom W.unop.hom.property.2.2 W.unop.hom.property.2.1 φ ψ
    (biratPullGamma G.core Z φ W) (biratPullGamma_preStep G.core Z φ W).2
    (biratPullGamma_preStep G.core Z φ W).1
    (biratPullAlpha G.core Z φ W) (biratPull_sq G.core Z φ W) _

/-! ## ★5. ★★★★関手 `𝒞^birat ⥤ 𝔽_{Φ^gp}` —— 原文の図式の中段 -/

/-- ★★★★**[FrdI] Proposition 4.4, (i)** —— **関手 `𝒞^birat → 𝔽_{Φ^gp}`**。

原文 (FrdI p.83):
> where the functors between elementary Frobenioids are those induced by the natural

★底と次数は下段(`biratToElem`)と同じ、★**零因子だけが新しい**。 -/
noncomputable def biratToElemGp : BiratCat P G ⥤ ElemFrobCat (Φ.gpOn (phi_integral P)) where
  obj A := ⟨(P.toElem.obj (biratDown P G A)).base⟩
  map f := ⟨biratBase f, biratDivGp f, biratDeg f⟩
  map_id A := by
    refine ElemFrobCat.Hom.ext ?_ ?_ ?_
    · show biratBase (toHomBirat (P := P) (G := G) (𝟙 (biratDown P G A))) = 𝟙 _
      rw [biratBase_toHomBirat, P.Base_id]
    · show biratDivGp (toHomBirat (P := P) (G := G) (𝟙 (biratDown P G A))) = 0
      rw [biratDivGp_toHomBirat, P.Div_id]
      exact toGp_zero
    · show biratDeg (toHomBirat (P := P) (G := G) (𝟙 (biratDown P G A))) = 1
      rw [biratDeg_toHomBirat, P.degFr_id]
  map_comp f g :=
    ElemFrobCat.Hom.ext (biratBase_comp G.core f g) (biratDivGp_comp f g)
      (biratDeg_comp G.core f g)

/-- ★★★**図式の上半分が可換**(対象) —— `𝒞 → 𝔽_Φ → 𝔽_{Φ^gp}` と
`𝒞 → 𝒞^birat → 𝔽_{Φ^gp}` は同じ対象を与える。 -/
theorem biratGp_square_obj (A : C) :
    (toBiratCat P G ⋙ biratToElemGp P G).obj A
      = (P.toElem ⋙ elemToGp P).obj A := rfl

variable {P G} in
/-- ★★★**図式の上半分が可換**(射) —— 底・零因子・次数のいずれも一致する。

★零因子が一致するのが `biratDivGp_toHomBirat`(`Div φ` の `Φ^gp` での像)である。 -/
theorem biratGp_square_map {A B : C} (φ : A ⟶ B) :
    (toBiratCat P G ⋙ biratToElemGp P G).map φ
      = (P.toElem ⋙ elemToGp P).map φ :=
  ElemFrobCat.Hom.ext (biratBase_toHomBirat φ) (biratDivGp_toHomBirat φ)
    (biratDeg_toHomBirat φ)

/-- ★locator —— **`Proposition 4.4, (i)`**。

★(i) の 3 つの主張はいずれも実装済である:
合成写像と圧 `𝒞^birat`(`BiratCat` / `biratCategory` / `compBirat`)、
自然な関手 `𝒞 → 𝒞^birat`(`toBiratCat`)、
そして **1-可換図式**(`biratToElemGp` と
`biratGp_square_obj` / `biratGp_square_map`)。 -/
def prop_4_4_i.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83, item := "Proposition 4.4, (i)",
    sectionId := "frdi-prop-4-4" }

/-! ## ★6. ★★★★`Φ^birat` —— `𝒪^×(A^birat)` の `Div^gp` による像

原文 (FrdI p.83):
> (iii) There exists a unique subfunctor of groups Φbirat ⊆Φgp such that

★★原文は `Φ^birat` を「中段が経由する最小の部分関手」として特徴づけ、
**各 `A^birat` で `𝒪^×(A^birat) ↠ Φ^birat(A^birat)` が全射**だと言う。
★したがって**像として定義すれば全射性は定義そのもの**になる。 -/

variable {P G} in
/-- ★`𝒞^birat` の合成の零因子(圏の `≫` の形)。 -/
theorem biratDivGp_comp' {A B E : BiratCat P G} (f : A ⟶ B) (g : B ⟶ E) :
    biratDivGp (f ≫ g) = gpMap _ (Φ.map (biratBase f)) (biratDivGp g)
      + ((biratDeg g : ℕ+) : ℕ) • biratDivGp f :=
  biratDivGp_comp f g

variable {P G} in
theorem biratDivGp_id (A : BiratCat P G) : biratDivGp (𝟙 A) = 0 := by
  show biratDivGp (toHomBirat (P := P) (G := G) (𝟙 (biratDown P G A))) = 0
  rw [biratDivGp_toHomBirat, P.Div_id]
  exact toGp_zero

variable {P G} in
/-- ★★**base-identity な射では `Φ^gp` の輸送が恒等になる**。 -/
theorem gpMap_biratBase_of_baseIdentity {A : BiratCat P G} {δ : A ⟶ A}
    (hδ : IsBaseIdentity (biratPre P G) δ)
    (x : Gp (Φ.val (P.toElem.obj (biratDown P G A)).base)) :
    gpMap _ (Φ.map (biratBase δ)) x = x := by
  have hb : biratBase δ = 𝟙 _ := by
    refine hδ.trans ?_
    show biratBase (toHomBirat (P := P) (G := G) (𝟙 (biratDown P G A))) = 𝟙 _
    rw [biratBase_toHomBirat, P.Base_id]
  rw [hb]
  have hmap : Φ.map (𝟙 ((P.toElem.obj (biratDown P G A)).base)) = AddMonoidHom.id _ := by
    ext y; exact Φ.map_id _ y
  rw [hmap, gpMap_id]
  rfl

variable {P G} in
/-- ★★**`𝒪^×` の上では `Div^gp` は加法的**。 -/
theorem biratDivGp_mul_otimes {A : BiratCat P G} {δ ε : End A}
    (hδ : δ ∈ OTimes (biratPre P G) A) (hε : ε ∈ OTimes (biratPre P G) A) :
    biratDivGp ((δ * ε : End A) : A ⟶ A)
      = biratDivGp (δ : A ⟶ A) + biratDivGp (ε : A ⟶ A) := by
  have h := biratDivGp_comp' (ε : A ⟶ A) (δ : A ⟶ A)
  rw [gpMap_biratBase_of_baseIdentity hε.1.1,
    show ((biratDeg (δ : A ⟶ A) : ℕ+) : ℕ) = 1 from by
      rw [show biratDeg (δ : A ⟶ A) = 1 from hδ.1.2]; rfl,
    one_smul] at h
  exact h

variable {P G} in
/-- ★**単元の逆元も `𝒪^×` に入る**。 -/
theorem otimes_inv_mem {A : BiratCat P G} {δ δ' : End A}
    (hδ : δ ∈ OTimes (biratPre P G) A)
    (h : (δ : A ⟶ A) ≫ (δ' : A ⟶ A) = 𝟙 A)
    (h' : (δ' : A ⟶ A) ≫ (δ : A ⟶ A) = 𝟙 A) : δ' ∈ OTimes (biratPre P G) A := by
  refine ⟨⟨?_, ?_⟩,
    (CategoryTheory.isUnit_iff_isIso (δ' : A ⟶ A)).mpr ⟨⟨(δ : A ⟶ A), h', h⟩⟩⟩
  · show (biratPre P G).Base (δ' : A ⟶ A) = (biratPre P G).Base (𝟙 A)
    have hb := congrArg (biratPre P G).Base h
    rw [(biratPre P G).Base_comp,
      show (biratPre P G).Base (δ : A ⟶ A) = (biratPre P G).Base (𝟙 A) from hδ.1.1,
      (biratPre P G).Base_id, Category.id_comp] at hb
    exact hb.trans ((biratPre P G).Base_id A).symm
  · show (biratPre P G).degFr (δ' : A ⟶ A) = 1
    have hd := congrArg (biratPre P G).degFr h
    rw [(biratPre P G).degFr_comp,
      show (biratPre P G).degFr (δ : A ⟶ A) = 1 from hδ.1.2,
      (biratPre P G).degFr_id, mul_one] at hd
    exact hd

/-- ★★★★**[FrdI] Proposition 4.4, (iii)** —— `Φ^birat(A_𝒟)`。

原文 (FrdI p.83):
> (iii) There exists a unique subfunctor of groups Φbirat ⊆Φgp such that

★**`𝒪^×(A^birat)` の `Div^gp` による像**として定義する。
★★全射性(原文の `𝒪^×(A^birat) ↠ Φ^birat(A^birat)`)は**定義そのもの**になる。 -/
noncomputable def phiBiratAt (A : BiratCat P G) :
    AddSubgroup (Gp (Φ.val (P.toElem.obj (biratDown P G A)).base)) where
  carrier := {x | ∃ δ ∈ OTimes (biratPre P G) A, biratDivGp (δ : A ⟶ A) = x}
  zero_mem' := ⟨1, (OTimes (biratPre P G) A).one_mem, biratDivGp_id A⟩
  add_mem' := by
    rintro x y ⟨δ, hδ, rfl⟩ ⟨ε, hε, rfl⟩
    exact ⟨δ * ε, (OTimes (biratPre P G) A).mul_mem hδ hε, biratDivGp_mul_otimes hδ hε⟩
  neg_mem' := by
    rintro x ⟨δ, hδ, rfl⟩
    obtain ⟨u, hu⟩ := hδ.2
    have h1 : (δ : A ⟶ A) ≫ ((↑u⁻¹ : End A) : A ⟶ A) = 𝟙 A := by
      have hv := u.inv_val
      rw [hu] at hv
      exact hv
    have h2 : ((↑u⁻¹ : End A) : A ⟶ A) ≫ (δ : A ⟶ A) = 𝟙 A := by
      have hv := u.val_inv
      rw [hu] at hv
      exact hv
    have hmem := otimes_inv_mem hδ h1 h2
    refine ⟨(↑u⁻¹ : End A), hmem, ?_⟩
    have hone : (δ * (↑u⁻¹ : End A) : End A) = 1 := h2
    have hsum := biratDivGp_mul_otimes hδ hmem
    rw [hone] at hsum
    have hz : biratDivGp ((1 : End A) : A ⟶ A) = 0 := biratDivGp_id A
    rw [hz] at hsum
    exact eq_neg_of_add_eq_zero_left ((add_comm _ _).trans hsum.symm)

/-! ## ★7. 残り —— `Φ^birat` の**関手性**

★★不変性(上)があるので、**余錐 `biratDivGpCocone` はそのまま組める**。
★その先の `map_comp`(合成則)が `𝔽_{Φ^gp}` への関手を与え、
`Φ^birat` はその像として定まる。
★★**`map_comp` は `Prop32.lean` の `pfDiv_comp` と同じ形の計算**になる見込みで、
そこでの測定(「目標の内側を書き換えようとしない」「`congrArg₂` でまとめて外側で `rw`」)が
そのまま効くはずである。
-/

end ABC3.Found.FrdI
