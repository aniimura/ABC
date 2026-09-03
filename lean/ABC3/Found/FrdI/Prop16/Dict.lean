/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop19
import ABC3.Found.FrdI.PlBkShuffle
import ABC3.Found.FrdI.Prop16.Restrict

/-!
# Prop16 —— 辞書——CfpCat の Frobenioid 構造の互換性

☆もとの 1 枚を**入れ子の切れ目**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2 u3 v3

variable {D : Type u} [Category.{v} D] {D' : Type u3} [Category.{v3} D']

section Dict

variable {C : Type u2} [Category.{v2} C] {Φ : MonoidOn.{v, u, w} D}

variable (P : PreFrobenioid C Φ) (G : D' ⥤ D)
  (hG : ∀ {A B : D'} (α : B ⟶ A), IsFSMMorphism α → IsFSMMorphism (G.map α))
  (hD' : IsTotallyEpimorphic D')
  (hcC : IsConnected (CfpCat P G)) (hcD' : IsConnected D')

theorem cfp_compat_Base {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD' hcC hcD').Base f = CfpCat.snd f := rfl

theorem cfp_compat_degFr {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD' hcC hcD').degFr f = P.degFr (CfpCat.fst f) := rfl

theorem cfp_compat_Div {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD' hcC hcD').Div f
      = Φ.map (@inv _ _ _ _ X.obj.hom X.property) (P.Div (CfpCat.fst f)) := rfl

/-- **(iii)** —— **Frobenius 次数**は射影で決まる。 -/
theorem cfp_degFr_eq {X Y : CfpCat P G} (f : X ⟶ Y) :
    (cfpPreFrobenioid P G hG hD' hcC hcD').degFr f = P.degFr (CfpCat.fst f) := rfl

/-- **(iv)** —— **base-isomorphism** は `𝒟'` 成分が同型であること。

★★**これが CFP の移送を支える一点である** —— `𝒟'` 成分の情報は
base-isomorphism の定義そのものに入っている。 -/
theorem cfp_baseIso_iff {X Y : CfpCat P G} (f : X ⟶ Y) :
    IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsIso (CfpCat.snd f) := Iff.rfl

/-- **(iii)** —— **linear** は射影で決まる。 -/
theorem cfp_linear_iff {X Y : CfpCat P G} (f : X ⟶ Y) :
    IsLinear (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsLinear P (CfpCat.fst f) := Iff.rfl

/-- **(iii)** —— **isometry** は射影で決まる。

★中身は「`Φ(α⁻¹)` が単射」の一点(`Definition 1.1, (ii), (a)`)。 -/
theorem cfp_isometric_iff {X Y : CfpCat P G} (f : X ⟶ Y) :
    IsIsometric (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsIsometric P (CfpCat.fst f) := by
  constructor
  · intro h
    exact Φ.map_injective (@inv _ _ _ _ X.obj.hom X.property) (h.trans (map_zero _).symm)
  · intro h
    show Φ.map (@inv _ _ _ _ X.obj.hom X.property) (P.Div (CfpCat.fst f)) = 0
    rw [show P.Div (CfpCat.fst f) = 0 from h]
    exact map_zero _

/-! ### ★**手順5**(2つの綴り問題への対処、確定版)

★`rw` を項で置き換えるだけでは足りない。**部分型の成分や `congrArg` の射影**は
綴りが食い違ったままなので、
**「綴りの決まった変数を `obtain ⟨h', rfl⟩ : ∃ h' : 正しい型, h' = h` で先に導入する」**
のが正しい対処である(`Prop19` の `Over.Hom.left` に使った定型と同じもの)。 -/

/-- ★`𝒞'` の射の `𝒞` 成分の底射は、`𝒟'` 成分と両端の同型で**完全に決まる**。

★★**これが CFP の移送を支える一点**である —— `𝒞` 成分の底射は自由ではない。 -/
theorem cfp_base_fst {X Y : CfpCat P G} (f : X ⟶ Y) [IsIso Y.obj.hom] :
    P.proj.map (CfpCat.fst f) = X.obj.hom ≫ G.map (CfpCat.snd f) ≫ inv Y.obj.hom := by
  rw [← Category.assoc, ← cfp_square f, Category.assoc, IsIso.hom_inv_id, Category.comp_id]

/-- **(iii)** —— `𝒞` の pull-back は `𝒞'` の pull-back(★**構成の向き**)。

★中身は「`𝒟'` 成分 `h` を与えると、`𝒞` 側で使うべき底射
`u = α_Z ≫ G(h) ≫ α_X⁻¹` が**一意に決まる**」の一点。 -/
theorem cfp_isPullBack_of {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsPullBack P (CfpCat.fst φ)) : IsPullBack (cfpPreFrobenioid P G hG hD' hcC hcD') φ := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  intro Z
  haveI hZ : IsIso Z.obj.hom := Z.property
  constructor
  · intro f₁ f₂ hf
    have hp := Subtype.ext_iff.mp hf
    have hs : CfpCat.snd f₁ = CfpCat.snd f₂ := congrArg Prod.snd hp
    have hcomp : (f₁ ≫ φ : Z ⟶ Y) = f₂ ≫ φ := congrArg Prod.fst hp
    have hc : CfpCat.fst f₁ ≫ CfpCat.fst φ = CfpCat.fst f₂ ≫ CfpCat.fst φ :=
      congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hcomp
    have hb : P.Base (CfpCat.fst f₁) = P.Base (CfpCat.fst f₂) := by
      show P.proj.map (CfpCat.fst f₁) = P.proj.map (CfpCat.fst f₂)
      rw [cfp_base_fst P G f₁, cfp_base_fst P G f₂, hs]
    exact InducedCategory.hom_ext
      (CommaMorphism.ext ((h Z.obj.left).1 (Subtype.ext (Prod.ext hc hb))) hs)
  · rintro ⟨⟨g, hh0⟩, hcond⟩
    -- ★手順5: 綴りの決まった変数を先に導入する
    obtain ⟨hh, rfl⟩ : ∃ hh : Z.obj.right ⟶ X.obj.right, hh = hh0 := ⟨hh0, rfl⟩
    have hcond' : CfpCat.snd g = hh ≫ CfpCat.snd φ := hcond
    obtain ⟨u, hu⟩ : ∃ u : P.proj.obj Z.obj.left ⟶ P.proj.obj X.obj.left,
        u = Z.obj.hom ≫ G.map hh ≫ inv X.obj.hom := ⟨_, rfl⟩
    have hbase : P.proj.map (CfpCat.fst g) = u ≫ P.proj.map (CfpCat.fst φ) := by
      rw [hu, cfp_base_fst P G g, cfp_base_fst P G φ, hcond', G.map_comp]
      simp only [Category.assoc]
      rw [← Category.assoc (inv X.obj.hom) X.obj.hom, IsIso.inv_hom_id, Category.id_comp]
    obtain ⟨f₁, hf₁⟩ := (h Z.obj.left).2 ⟨(CfpCat.fst g, u), hbase⟩
    have hp := Subtype.ext_iff.mp hf₁
    have h1 : (f₁ ≫ CfpCat.fst φ) = CfpCat.fst g := congrArg Prod.fst hp
    have h2 : P.proj.map f₁ = u := congrArg Prod.snd hp
    have hw : P.proj.map f₁ ≫ X.obj.hom = Z.obj.hom ≫ G.map hh := by
      rw [h2, hu]
      simp only [Category.assoc]
      rw [IsIso.inv_hom_id, Category.comp_id]
    refine ⟨InducedCategory.homMk ⟨f₁, hh, hw⟩, Subtype.ext (Prod.ext ?_ rfl)⟩
    exact InducedCategory.hom_ext (CommaMorphism.ext h1 hcond'.symm)

/-- ★`𝒞'` の射は、両成分が同型なら同型。 -/
theorem cfp_isIso_of {X Y : CfpCat P G} (f : X ⟶ Y) (h1 : IsIso (CfpCat.fst f))
    (h2 : IsIso (CfpCat.snd f)) : IsIso f := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  haveI := h1
  haveI := h2
  have hw : P.proj.map (inv (CfpCat.fst f)) ≫ X.obj.hom
      = Y.obj.hom ≫ G.map (inv (CfpCat.snd f)) := by
    rw [P.proj.map_inv, G.map_inv, IsIso.inv_comp_eq, ← Category.assoc, IsIso.eq_comp_inv]
    exact (cfp_square f).symm
  refine ⟨InducedCategory.homMk ⟨inv (CfpCat.fst f), inv (CfpCat.snd f), hw⟩, ?_, ?_⟩
  · exact InducedCategory.hom_ext (CommaMorphism.ext (IsIso.hom_inv_id _) (IsIso.hom_inv_id _))
  · exact InducedCategory.hom_ext (CommaMorphism.ext (IsIso.inv_hom_id _) (IsIso.inv_hom_id _))

/-- ★`𝒞'` の base-isomorphism は `𝒞` の base-isomorphism。

★`cfp_base_fst` により `𝒞` 成分の底射は同型3つの合成になる。 -/
theorem cfp_baseIso_fst {X Y : CfpCat P G} (f : X ⟶ Y)
    (h : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f) :
    IsBaseIsomorphism P (CfpCat.fst f) := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  haveI : IsIso (CfpCat.snd f) := h
  show IsIso (P.proj.map (CfpCat.fst f))
  rw [cfp_base_fst P G f]
  infer_instance

/-- **(iii)** —— **co-angular** は `𝒞` から `𝒞'` へ降りる(★構成の向き)。 -/
theorem cfp_coAngular_of {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsCoAngular P (CfpCat.fst φ)) :
    IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ := by
  intro Z W γ β α hfac hαl hβi hβs hdisj
  have hfac' : CfpCat.fst φ = CfpCat.fst γ ≫ CfpCat.fst β ≫ CfpCat.fst α :=
    congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hfac
  have hβsC : IsPreStep P (CfpCat.fst β) :=
    ⟨hβs.1, cfp_baseIso_fst P G hG hD' hcC hcD' β hβs.2⟩
  have hdisjC : IsBaseIsomorphism P (CfpCat.fst α) ∨ IsBaseIsomorphism P (CfpCat.fst γ) :=
    hdisj.imp (cfp_baseIso_fst P G hG hD' hcC hcD' α) (cfp_baseIso_fst P G hG hD' hcC hcD' γ)
  have := h _ _ (CfpCat.fst γ) (CfpCat.fst β) (CfpCat.fst α) hfac' hαl
    ((cfp_isometric_iff P G hG hD' hcC hcD' β).mp hβi) hβsC hdisjC
  exact cfp_isIso_of P G β this hβs.2

/-- ★`𝒞` の対象を、底が `G` の像と同型であるときに `𝒞'` へ持ち上げる。 -/
def cfpMk (A : C) (W : D') (e : P.proj.obj A ⟶ G.obj W) [IsIso e] : CfpCat P G :=
  ⟨⟨A, W, e⟩, inferInstanceAs (IsIso e)⟩

/-- ★`𝒞'` の射を両成分と四角形から作る。 -/
def cfpHom {X Y : CfpCat P G} (u : X.obj.left ⟶ Y.obj.left) (v : X.obj.right ⟶ Y.obj.right)
    (w : P.proj.map u ≫ Y.obj.hom = X.obj.hom ≫ G.map v) : X ⟶ Y :=
  InducedCategory.homMk ⟨u, v, w⟩

/-- ★`𝒞'` の同型の `𝒞` 成分は同型。 -/
theorem cfp_isIso_fst {X Y : CfpCat P G} (f : X ⟶ Y) (h : IsIso f) : IsIso (CfpCat.fst f) := by
  obtain ⟨g, h1, h2⟩ := h.out
  refine ⟨CfpCat.fst g, ?_, ?_⟩
  · exact congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) h1
  · exact congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) h2

/-! ### ★**手順6**(2つの綴り問題への対処、3つ目)

★**射影が簡約されない `def` を作らず、構造リテラルを直接書く**。
`cfpMk` のような補助定義を挟むと `(cfpMk …).obj.right` が `Y.obj.right` に
**簡約されず**、`rw` も インスタンス合成も落ちる。
★手順5(部分型の成分)では届かない、3つ目の出方である。 -/

/-- **(iii)** —— **co-angular** は `𝒞'` から `𝒞` へ上がる(★**難しい向き**)。

★★**要点**: `co-angular` の定義に入っている**選言**
「`α` か `γ` が base-isomorphism」が、**中間対象を `𝒞'` へ持ち上げる橋**になる。
`β` は pre-step なので必ず base-isomorphism であり、
選言の側からもう一方の端まで **base-isomorphism の鎖**が繋がるので、
`Base Z₀` も `Base W₀` も `G` の像と同型になる。

★**`G` の本質的全射性は要らない** —— 使うのは「与えられた分解の中の base-isomorphism」だけ。
★★**これが CFP 版の仕分け基準**である: **定義の中に base-isomorphism の鎖があるか**。 -/
theorem cfp_coAngular_to {X Y : CfpCat P G} (φ : X ⟶ Y)
    (h : IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ) : IsCoAngular P (CfpCat.fst φ) := by
  haveI hX : IsIso X.obj.hom := X.property
  haveI hY : IsIso Y.obj.hom := Y.property
  intro Z₀ W₀ γ₀ β₀ α₀ hfac hα₀l hβ₀i hβ₀s hdisj
  haveI hβb : IsIso (P.proj.map β₀) := hβ₀s.2
  have hbf : P.proj.map (CfpCat.fst φ) ≫ Y.obj.hom = X.obj.hom ≫ G.map (CfpCat.snd φ) :=
    cfp_square φ
  have hfacp : P.proj.map (CfpCat.fst φ)
      = P.proj.map γ₀ ≫ P.proj.map β₀ ≫ P.proj.map α₀ := by
    rw [hfac, P.proj.map_comp, P.proj.map_comp]
  rcases hdisj with hα | hγ
  · -- ★`α₀` が base-isomorphism: `Y` の側から鎖を辿る
    haveI hαb : IsIso (P.proj.map α₀) := hα
    have hz : IsIso (P.proj.map β₀ ≫ P.proj.map α₀ ≫ Y.obj.hom) := inferInstance
    have hw : IsIso (P.proj.map α₀ ≫ Y.obj.hom) := inferInstance
    refine cfp_isIso_fst P G
      (X := ⟨⟨Z₀, Y.obj.right, P.proj.map β₀ ≫ P.proj.map α₀ ≫ Y.obj.hom⟩, hz⟩)
      (Y := ⟨⟨W₀, Y.obj.right, P.proj.map α₀ ≫ Y.obj.hom⟩, hw⟩)
      (InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩) ?_
    refine h ⟨⟨Z₀, Y.obj.right, P.proj.map β₀ ≫ P.proj.map α₀ ≫ Y.obj.hom⟩, hz⟩
      ⟨⟨W₀, Y.obj.right, P.proj.map α₀ ≫ Y.obj.hom⟩, hw⟩
      (InducedCategory.homMk ⟨γ₀, (CfpCat.snd φ : X.obj.right ⟶ Y.obj.right), ?_⟩)
      (InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩)
      (InducedCategory.homMk ⟨α₀, 𝟙 _, by simp⟩) ?_
      hα₀l ((cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hβ₀i) ⟨hβ₀s.1, by show IsIso (𝟙 _); infer_instance⟩
      (Or.inl (by show IsIso (𝟙 _); infer_instance))
    · show P.proj.map γ₀ ≫ P.proj.map β₀ ≫ P.proj.map α₀ ≫ Y.obj.hom
        = X.obj.hom ≫ G.map (CfpCat.snd φ)
      rw [← hbf, hfacp]
      simp only [Category.assoc]
    · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
      show CfpCat.snd φ = CfpCat.snd φ ≫ 𝟙 _ ≫ 𝟙 _
      simp
  · -- ★`γ₀` が base-isomorphism: `X` の側から鎖を辿る
    haveI hγb : IsIso (P.proj.map γ₀) := hγ
    have hz : IsIso (inv (P.proj.map γ₀) ≫ X.obj.hom) := inferInstance
    have hw : IsIso (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ X.obj.hom) := inferInstance
    refine cfp_isIso_fst P G
      (X := ⟨⟨Z₀, X.obj.right, inv (P.proj.map γ₀) ≫ X.obj.hom⟩, hz⟩)
      (Y := ⟨⟨W₀, X.obj.right,
        inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ X.obj.hom⟩, hw⟩)
      (InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩) ?_
    refine h ⟨⟨Z₀, X.obj.right, inv (P.proj.map γ₀) ≫ X.obj.hom⟩, hz⟩
      ⟨⟨W₀, X.obj.right,
        inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ X.obj.hom⟩, hw⟩
      (InducedCategory.homMk ⟨γ₀, 𝟙 _, by simp⟩)
      (InducedCategory.homMk ⟨β₀, 𝟙 _, by simp⟩)
      (InducedCategory.homMk ⟨α₀, (CfpCat.snd φ : X.obj.right ⟶ Y.obj.right), ?_⟩) ?_
      hα₀l ((cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hβ₀i) ⟨hβ₀s.1, by show IsIso (𝟙 _); infer_instance⟩
      (Or.inr (by show IsIso (𝟙 _); infer_instance))
    · show P.proj.map α₀ ≫ Y.obj.hom
        = (inv (P.proj.map β₀) ≫ inv (P.proj.map γ₀) ≫ X.obj.hom) ≫ G.map (CfpCat.snd φ)
      rw [Category.assoc, Category.assoc, ← hbf, hfacp]
      simp only [Category.assoc]
      rw [← Category.assoc (inv (P.proj.map γ₀)), IsIso.inv_hom_id, Category.id_comp,
        ← Category.assoc (inv (P.proj.map β₀)), IsIso.inv_hom_id, Category.id_comp]
    · refine InducedCategory.hom_ext (CommaMorphism.ext hfac ?_)
      show CfpCat.snd φ = 𝟙 _ ≫ 𝟙 _ ≫ CfpCat.snd φ
      simp

/-- **(iii)** —— co-angular は射影で決まる(両向き)。 -/
theorem cfp_coAngular_iff {X Y : CfpCat P G} (φ : X ⟶ Y) :
    IsCoAngular (cfpPreFrobenioid P G hG hD' hcC hcD') φ ↔ IsCoAngular P (CfpCat.fst φ) :=
  ⟨cfp_coAngular_to P G hG hD' hcC hcD' φ, cfp_coAngular_of P G hG hD' hcC hcD' φ⟩

/-- **(iii)** —— **LB-invertible** は射影で決まる。

★`LB-invertible = co-angular ∧ isometric` なので、上の2つの合成。 -/
theorem cfp_lbInvertible_iff {X Y : CfpCat P G} (φ : X ⟶ Y) :
    IsLBInvertible (cfpPreFrobenioid P G hG hD' hcC hcD') φ ↔ IsLBInvertible P (CfpCat.fst φ) :=
  and_congr (cfp_coAngular_iff P G hG hD' hcC hcD' φ) (cfp_isometric_iff P G hG hD' hcC hcD' φ)

/-! ### ★★**分類表** —— 「型は等しいが字面が違う」問題の症状・原因・対処

★**個数ではなく「症状が分岐したこと」が表を作る合図**だった。
同じ `def` 由来でも対処が違うと分かった時点で分ける。
★**切り替え点(新しいファイル/section/新しい種類の補題)でこの表を引く。**

| # | 症状(出るエラー) | 原因 | 対処 |
|---|---|---|---|
| 1 | `rw` が「目標に見えているのに見つからない」 | 同じ型に2つの綴り(`def` が展開されない) | **`rw` を使わず `Eq.trans`/`congrArg` で項として繋ぐ** |
| 2 | `inv` を書くと `motive is not type correct` | `InducedCategory.Hom` に包まれて型が簡約されない | **`inv` を書かず `IsIso.out` から逆射を取る**(例外: 主張自体に `inv` が現れるときは `IsIso.eq_inv_comp` にインスタンスを明示) |
| 3 | 部分型 `{p : _ × _ // _}` の成分の型が食い違う | 成分の型が別綴りで書かれている | **`obtain ⟨h', rfl⟩ : ∃ h' : 正しい型, h' = h` で綴りの決まった変数を先に導入** |
| 4 | `(f …).obj.right` が簡約されず `rw`/合成が落ちる | 補助 `def` の射影は簡約されない | **補助 `def` を作らず構造リテラルを直接書く** |
| 5 | `failed to synthesize instance` (文脈に `haveI` があるのに) | インスタンス型が `def`(例: `cfpProp`)で書かれ展開されない | **`have h : IsIso … := inferInstance` を先に置き、その項を明示的に渡す** |
| 6 | `refine h _ _ …` の穴が埋まらない | 中間対象がメタ変数のまま残る | **対象を明示的に書く**(`_` に頼らない) |
| 7 | `include` し忘れで `Unknown identifier` | 主張に現れない変数は自動包含されない | **最初に `include F in`** |
-/

/-! ## ★(iv) —— base-isomorphism についての3クラスと `𝒪^▷`

原文 (FrdI p.28):
> tively, pre-step; step) if and only if its projection to C is. Moreover, the projection

★★**基準の4例目**: (iv) が「**base-isomorphism について**」と限定して述べられているのは、
**その仮定が鎖そのものを与える**からである。仮定を外すと
「`fst f` が base-iso でも `snd f` が同型とは限らない」(`G` は同型を反映しない)ので、
**⟸ 向きが壊れる**。★**原文の限定は必要であり、我々の基準がその理由を説明する。** -/

/-- **(iv)** —— base-isomorphism については **pre-step** が射影で決まる。 -/
theorem cfp_preStep_iff {X Y : CfpCat P G} (f : X ⟶ Y)
    (hb : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f) :
    IsPreStep (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsPreStep P (CfpCat.fst f) :=
  ⟨fun h => ⟨h.1, cfp_baseIso_fst P G hG hD' hcC hcD' f hb⟩, fun h => ⟨h.1, hb⟩⟩

/-- **(iv)** —— base-isomorphism については **Frobenius 型**が射影で決まる。 -/
theorem cfp_frobType_iff {X Y : CfpCat P G} (f : X ⟶ Y)
    (hb : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f) :
    IsFrobeniusType (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsFrobeniusType P (CfpCat.fst f) :=
  ⟨fun h => ⟨(cfp_lbInvertible_iff P G hG hD' hcC hcD' f).mp h.1, cfp_baseIso_fst P G hG hD' hcC hcD' f hb⟩,
   fun h => ⟨(cfp_lbInvertible_iff P G hG hD' hcC hcD' f).mpr h.1, hb⟩⟩

/-- **(iv)** —— base-isomorphism については **step** が射影で決まる。

★同型性の両向きは `cfp_isIso_of` / `cfp_isIso_fst`。 -/
theorem cfp_step_iff {X Y : CfpCat P G} (f : X ⟶ Y)
    (hb : IsBaseIsomorphism (cfpPreFrobenioid P G hG hD' hcC hcD') f) :
    IsStep (cfpPreFrobenioid P G hG hD' hcC hcD') f ↔ IsStep P (CfpCat.fst f) := by
  constructor
  · rintro ⟨hs, hni⟩
    exact ⟨(cfp_preStep_iff P G hG hD' hcC hcD' f hb).mp hs,
      fun hi => hni (cfp_isIso_of P G f hi hb)⟩
  · rintro ⟨hs, hni⟩
    exact ⟨(cfp_preStep_iff P G hG hD' hcC hcD' f hb).mpr hs,
      fun hi => hni (cfp_isIso_fst P G f hi)⟩

/-- ★`𝒞'` の base-identity 自己射は、`𝒟'` 成分が `𝟙` であること。 -/
theorem cfp_baseIdentity_iff {A : CfpCat P G} (e : A ⟶ A) :
    IsBaseIdentity (cfpPreFrobenioid P G hG hD' hcC hcD') e ↔ CfpCat.snd e = 𝟙 A.obj.right :=
  Iff.rfl

/-- ★その `𝒞` 成分は `𝒞` の base-identity。 -/
theorem cfp_baseIdentity_fst {A : CfpCat P G} (e : A ⟶ A)
    (h : IsBaseIdentity (cfpPreFrobenioid P G hG hD' hcC hcD') e) :
    IsBaseIdentity P (CfpCat.fst e) := by
  haveI hA : IsIso A.obj.hom := A.property
  have hsq := cfp_square e
  rw [(cfp_baseIdentity_iff P G hG hD' hcC hcD' e).mp h, G.map_id, Category.comp_id] at hsq
  show P.proj.map (CfpCat.fst e) = P.proj.map (𝟙 A.obj.left)
  rw [P.proj.map_id]
  exact (cancel_mono A.obj.hom).mp (hsq.trans (Category.id_comp _).symm)

/-- **(iv)** —— 射影は **`𝒪^▷` のモノイド同型**を誘導する。

★原文の当該行(`functor C′ →C determines a bijection of monoids O▷(A′)`)は
**`′`(prime) と `▷` を含むため逐語照合できない**。書き換えず、照合できない事実として記す。

★★**基準の5例目**: `base-identity` の定義が **`𝒟'` 成分を `𝟙` に固定する**ので、
`𝒪^▷` の元は `𝒞` 成分だけで決まる。**鎖どころか、`𝒟'` 成分が一意に定まる。** -/
def cfpOTriEquiv (A : CfpCat P G) :
    OTri (cfpPreFrobenioid P G hG hD' hcC hcD') A ≃* OTri P A.obj.left where
  toFun e := ⟨CfpCat.fst (e : End A), cfp_baseIdentity_fst P G hG hD' hcC hcD' _ e.2.1, e.2.2⟩
  invFun e :=
    ⟨InducedCategory.homMk ⟨(e : End A.obj.left), 𝟙 _, by
      have h1 : P.proj.map (e : End A.obj.left) = 𝟙 _ := e.2.1.trans (P.Base_id _)
      rw [h1, G.map_id, Category.comp_id, Category.id_comp]⟩,
     by show CfpCat.snd _ = 𝟙 _; rfl, e.2.2⟩
  left_inv e := by
    refine Subtype.ext (InducedCategory.hom_ext (CommaMorphism.ext rfl ?_))
    exact ((cfp_baseIdentity_iff P G hG hD' hcC hcD' (e : End A)).mp e.2.1).symm
  right_inv e := Subtype.ext rfl
  map_mul' x y := rfl

/-! ## ★(v) —— 対象タイプ

原文 (FrdI p.28):
> and only if it projects to such an object of C.

★★**仕分け**(基準を対象タイプ向けに変形したもの):

* **(a) 鎖型** —— 定義に現れる射が base-isomorphism なので中間対象が持ち上がる。
  `isotropic`(isometric pre-step) / `sub-quasi-Frobenius-trivial`(pre-step) /
  `Frobenius-isotropic`(Frobenius 型) / `perfect`(仮定の `base-isomorphic` が鎖)
* **(b) `𝒟'` 成分固定型** —— `base-identity` が `𝒟'` 成分を `𝟙` に固定する。
  `Frobenius-trivial` / `quasi-Frobenius-trivial` / `Frobenius-normalized` /
  `unit-trivial`(`𝒪^×` 経由) / `group-like`(`Φ` の同型経由)
* ★**(c) 未解決** —— `metrically trivial` / `base-trivial`。
  **結論が「`Nonempty (X ≅ A)`」**で、`𝒞'` の同型を作るには四角形が
  `𝒞` 成分の `Base` を1つに指定してしまう。**`Aut-ample` 相当**が要り、仮定にない。
-/

/-- **(v)** —— **isotropic** は射影で決まる(★(a) 鎖型)。 -/
theorem cfp_isotropic_iff (A : CfpCat P G) :
    IsIsotropic (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔ IsIsotropic P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · intro h Dd₀ φ₀ hi hs
    haveI hφb : IsIso (P.proj.map φ₀) := hs.2
    have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
    refine cfp_isIso_fst P G
      (X := A) (Y := ⟨⟨Dd₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩)
      (InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩) ?_
    refine h ⟨⟨Dd₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩
      (InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩)
      ((cfp_isometric_iff P G hG hD' hcC hcD' _).mpr hi) ⟨hs.1, by show IsIso (𝟙 _); infer_instance⟩
  · intro h Dd' φ hi hs
    exact cfp_isIso_of P G φ
      (h Dd'.obj.left (CfpCat.fst φ) ((cfp_isometric_iff P G hG hD' hcC hcD' φ).mp hi)
        ((cfp_preStep_iff P G hG hD' hcC hcD' φ hs.2).mp hs)) hs.2

/-! ### ★★(c) の 2 件 —— **片向きだけ実装する**(2026-08-16)

★`metrically trivial` / `base-trivial` は結論が `Nonempty (X ≅ A)` である。

★★**`⟸`(`𝒞` で成り立つ ⟹ `𝒞'` で成り立つ)は出ない** ——
`𝒞'` の同型を作るには四角形
`A.hom ∘ Base f = G(g) ∘ Dd'.hom` を満たす `(f, g)` が要り、
★**`f` の `Base` が 1 つに指定されてしまう。**
`base-trivial` が与えるのは**ある**同型であって、底を指定した同型ではない。
★**直すには `Aut-ample` 相当(または `G` の充満性)が要るが、`Proposition 1.6` の
仮定は「FSM 射を FSM 射に送る関手」だけである。**

★★**`⟹` は出る。** `𝒞` の対象 `Dd₀` を `A.obj.right` と組んで `𝒞'` の対象にし、
`𝒞'` の主張を当てて `fst` を取ればよい ——
★**この向きは「底を指定する」必要がない。**

★★**原文が (vi) を「if」(片向き)としか書いていないことに注意** ——
★**著者は (v) と (vi) で向きを書き分けている。**
★**(v) の 2 件については、我々は片向きしか出せていない。**
-/

/-- ★**(v) の `metrically trivial`(★片向き)**。 -/
theorem cfp_metricallyTrivial_mp (A : CfpCat P G)
    (h : IsMetricallyTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A) :
    IsMetricallyTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  intro Dd₀ φ₀ hc hs
  haveI hφb : IsIso (P.proj.map φ₀) := hs.2
  have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
  obtain ⟨w⟩ := h ⟨⟨Dd₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩
    (InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩)
    ((cfp_coAngular_iff P G hG hD' hcC hcD' _).mpr hc)
    ⟨hs.1, by show IsIso (𝟙 _); infer_instance⟩
  haveI := cfp_isIso_fst P G w.hom inferInstance
  exact ⟨asIso (CfpCat.fst w.hom)⟩

/-- ★**(v) の `base-trivial`(★片向き)**。 -/
theorem cfp_baseTrivial_mp (A : CfpCat P G)
    (h : IsBaseTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A) :
    IsBaseTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  intro Dd₀ hbi
  obtain ⟨e⟩ := hbi
  have hh : (P.toElem.obj Dd₀).base ≅ G.obj A.obj.right :=
    e.symm ≪≫ (@asIso _ _ _ _ A.obj.hom hA)
  obtain ⟨w⟩ := h ⟨⟨Dd₀, A.obj.right, hh.hom⟩, inferInstanceAs (IsIso hh.hom)⟩
    ⟨Iso.refl _⟩
  haveI := cfp_isIso_fst P G w.hom inferInstance
  exact ⟨asIso (CfpCat.fst w.hom)⟩

/-! ### ★★(c) の `⟸` —— **`Aut-ample` を足せば通る**(2026-08-16)

★**両向きが 1 点に帰着することを先に確定させた**:
`𝒞'` の同型は対 `(θ, g)` で四角形
`P.proj.map θ ≫ A.obj.hom = X.obj.hom ≫ G.map g` を満たすものなので、
★**`θ` の底が 1 つに指定される。**

* `base-trivial`: 指定される底は `X.obj.hom ≫ G.map g ≫ inv A.obj.hom`
* `metrically trivial`: `cfp_square` を解くと `inv (P.Base (fst φ))`

★★**したがって要るのはただ 1 つ**:
> **(†) 指定された同型 `v : Base X₀ ≅ Base A₀` に対し、`Base θ = v` なる同型 `θ : X₀ ≅ A₀`**

★`X₀ := A₀` と取れるので **(†) ⟺ `A₀` が `Aut-ample`** である。
★**`Definition 1.3` の 21 条からは (†) は出ない**(測定の記録は上の docstring)。
★★**ここでは原文の仮定に `Aut-ample` を足して閉じる。** -/

/-- ★★**底を指定した同型**(`Aut-ample` の言い換え) ——
「**ある**同型」を「**指定の底を持つ**同型」に直す。

★`base-trivial` / `metrically trivial` はどちらも「ある同型」しか与えないので、
★**この 1 本が (v) の残り 2 件の共通の心臓部**である。 -/
theorem cfp_iso_of_isAutAmple {A₀ X₀ : C} (haa : IsAutAmple P A₀)
    (θ₀ : X₀ ≅ A₀) (v : P.proj.obj X₀ ⟶ P.proj.obj A₀) [IsIso v] :
    ∃ θ : X₀ ⟶ A₀, IsIso θ ∧ P.proj.map θ = v := by
  haveI hb0 : IsIso (P.proj.map θ₀.hom) := by
    refine ⟨P.proj.map θ₀.inv, ?_, ?_⟩
    · rw [← P.proj.map_comp, θ₀.hom_inv_id, P.proj.map_id]
    · rw [← P.proj.map_comp, θ₀.inv_hom_id, P.proj.map_id]
  have hviso : IsIso (P.proj.map θ₀.inv ≫ v) := by
    haveI : IsIso (P.proj.map θ₀.inv) := by
      refine ⟨P.proj.map θ₀.hom, ?_, ?_⟩
      · rw [← P.proj.map_comp, θ₀.inv_hom_id, P.proj.map_id]
      · rw [← P.proj.map_comp, θ₀.hom_inv_id, P.proj.map_id]
    infer_instance
  obtain ⟨c, hci, hcb⟩ := haa (P.proj.map θ₀.inv ≫ v) hviso
  haveI := hci
  refine ⟨θ₀.hom ≫ (c : A₀ ⟶ A₀), inferInstance, ?_⟩
  rw [P.proj.map_comp,
    show P.proj.map (c : A₀ ⟶ A₀) = P.proj.map θ₀.inv ≫ v from hcb,
    ← Category.assoc, ← P.proj.map_comp, θ₀.hom_inv_id, P.proj.map_id, Category.id_comp]

/-- ★★**(v) の `metrically trivial`(`⟸`、`Aut-ample` 仮定)**。

★`𝒞'` の co-angular pre-step `φ` の `𝒞` 成分に metric triviality を当てて同型 `κ` を得、
★**`cfp_iso_of_isAutAmple` で底を `inv (Base (fst φ))` に直す。**
`𝒟'` 成分は `inv (snd φ)` を取れば四角形が `cfp_square` から出る。 -/
theorem cfp_metricallyTrivial_mpr (A : CfpCat P G) (haa : IsAutAmple P A.obj.left)
    (h : IsMetricallyTrivial P A.obj.left) :
    IsMetricallyTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
  haveI hA : IsIso A.obj.hom := A.property
  intro Dd' φ hc hs
  haveI hD : IsIso Dd'.obj.hom := Dd'.property
  have hψc : IsCoAngular P (CfpCat.fst φ) := (cfp_coAngular_iff P G hG hD' hcC hcD' φ).mp hc
  have hψs : IsPreStep P (CfpCat.fst φ) := (cfp_preStep_iff P G hG hD' hcC hcD' φ hs.2).mp hs
  haveI hψb : IsIso (P.proj.map (CfpCat.fst φ)) := hψs.2
  haveI hsb : IsIso (CfpCat.snd φ) := hs.2
  obtain ⟨κ⟩ := h Dd'.obj.left (CfpCat.fst φ) hψc hψs
  obtain ⟨θ, hθi, hθv⟩ :=
    cfp_iso_of_isAutAmple P haa κ (inv (P.proj.map (CfpCat.fst φ)))
  haveI := hθi
  have hsq : P.proj.map θ ≫ A.obj.hom
      = Dd'.obj.hom ≫ G.map (inv (CfpCat.snd φ)) := by
    rw [hθv, G.map_inv, IsIso.eq_comp_inv, Category.assoc, ← cfp_square φ,
      ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
  haveI : IsIso (cfpHom P G (X := Dd') (Y := A) θ (inv (CfpCat.snd φ)) hsq) :=
    cfp_isIso_of P G _ hθi (inferInstanceAs (IsIso (inv (CfpCat.snd φ))))
  exact ⟨asIso (cfpHom P G (X := Dd') (Y := A) θ (inv (CfpCat.snd φ)) hsq)⟩

/-- ★★**(v) の `base-trivial`(`⟸`、`Aut-ample` 仮定)**。

★`𝒟'` 成分の同型 `e'` を `G` で送って `𝒞` 成分の底同型 `v` を作り、
base-triviality で**ある**同型を取ってから
★**`cfp_iso_of_isAutAmple` で底を `v` に直す。** -/
theorem cfp_baseTrivial_mpr (A : CfpCat P G) (haa : IsAutAmple P A.obj.left)
    (h : IsBaseTrivial P A.obj.left) :
    IsBaseTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
  haveI hA : IsIso A.obj.hom := A.property
  intro X hbi
  haveI hX : IsIso X.obj.hom := X.property
  obtain ⟨e'⟩ := hbi
  have he' : A.obj.right ≅ X.obj.right := e'
  haveI hv : IsIso (X.obj.hom ≫ G.map he'.inv ≫ inv A.obj.hom) := inferInstance
  obtain ⟨θ₀⟩ := h X.obj.left
    ⟨(@asIso _ _ _ _ (X.obj.hom ≫ G.map he'.inv ≫ inv A.obj.hom) hv).symm⟩
  obtain ⟨θ, hθi, hθv⟩ :=
    cfp_iso_of_isAutAmple P haa θ₀ (X.obj.hom ≫ G.map he'.inv ≫ inv A.obj.hom)
  haveI := hθi
  have hsq : P.proj.map θ ≫ A.obj.hom = X.obj.hom ≫ G.map he'.inv := by
    rw [hθv, Category.assoc, Category.assoc, IsIso.inv_hom_id, Category.comp_id]
  haveI : IsIso (cfpHom P G (X := X) (Y := A) θ he'.inv hsq) :=
    cfp_isIso_of P G _ hθi (inferInstanceAs (IsIso he'.inv))
  exact ⟨asIso (cfpHom P G (X := X) (Y := A) θ he'.inv hsq)⟩

/-- ★★★**上の 2 本の仮定は空虚ではない**。

★★**「仮定を足せば何でも証明できる」への歯止め**である ——
足した `Aut-ample` が `base-trivial` と両立しないなら上の 2 本は**空虚に真**になるが、
★**両方を満たす対象は実在する**: `𝔽_Φ` の対象はすべてそうである
(`Proposition 1.5, (i)`、`elemFrob_autAmple` / `elemFrob_baseTrivial`)。

★あわせて記録: 上の 2 本は `haa` と `h` を**どちらも使う**(片方だけでは通らない)。
★そして `.src` は**付けない** —— `Proposition 1.6` の (v) を完全に実装したという
主張ではないからである。穴は `Gap/FrdI/Section1.lean` の `Gap_1_6_v` のまま。 -/
theorem autAmple_and_baseTrivial_nonvacuous {Φ₀ : MonoidOn.{v, u, w} D}
    (hD₀ : IsTotallyEpimorphic D) [IsConnected D]
    (hpd₀ : ∀ A : D, IsPreDivisorial (Φ₀.val A)) (A : ElemFrobCat Φ₀) :
    IsAutAmple (elemPreFrobenioid Φ₀ hD₀ hpd₀) A ∧
      IsBaseTrivial (elemPreFrobenioid Φ₀ hD₀ hpd₀) A :=
  ⟨elemFrob_autAmple Φ₀ hD₀ hpd₀ A, elemFrob_baseTrivial Φ₀ hD₀ hpd₀ A⟩

/-- **(v)** —— **Frobenius-isotropic** は射影で決まる(★(a) 鎖型)。 -/
theorem cfp_frobIsotropic_iff (A : CfpCat P G) :
    IsFrobeniusIsotropic (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔
      IsFrobeniusIsotropic P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · rintro ⟨Dd', φ, hft, hiso⟩
    exact ⟨Dd'.obj.left, CfpCat.fst φ, (cfp_frobType_iff P G hG hD' hcC hcD' φ hft.2).mp hft,
      (cfp_isotropic_iff P G hG hD' hcC hcD' Dd').mp hiso⟩
  · rintro ⟨Dd₀, φ₀, hft, hiso⟩
    haveI hφb : IsIso (P.proj.map φ₀) := hft.2
    have hzi : IsIso (inv (P.proj.map φ₀) ≫ A.obj.hom) := inferInstance
    refine ⟨⟨⟨Dd₀, A.obj.right, inv (P.proj.map φ₀) ≫ A.obj.hom⟩, hzi⟩,
      InducedCategory.homMk ⟨φ₀, 𝟙 _, by simp⟩, ?_, ?_⟩
    · exact (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr hft
    · exact (cfp_isotropic_iff P G hG hD' hcC hcD' _).mpr hiso

/-- **(v)** —— **quasi-Frobenius-trivial** は射影で決まる(★(b) `𝒟'` 成分固定型)。 -/
theorem cfp_quasiFrobTrivial_iff (A : CfpCat P G) :
    IsQuasiFrobeniusTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔
      IsQuasiFrobeniusTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · intro h n
    obtain ⟨φ, hbi, hdeg⟩ := h n
    exact ⟨CfpCat.fst φ, cfp_baseIdentity_fst P G hG hD' hcC hcD' φ hbi, hdeg⟩
  · intro h n
    obtain ⟨φ₀, hbi, hdeg⟩ := h n
    have hid : P.proj.map φ₀ = 𝟙 _ := hbi.trans (P.Base_id _)
    exact ⟨InducedCategory.homMk ⟨φ₀, 𝟙 _, by rw [hid, G.map_id, Category.comp_id,
      Category.id_comp]⟩, by show CfpCat.snd _ = 𝟙 _; rfl, hdeg⟩

/-- **(v)** —— **Frobenius-trivial** は射影で決まる(★(b) 型)。

★`ζ : ℕ≥1 →* End A` を運ぶとき、**`𝒟'` 成分をすべて `𝟙` に取る**ので
モノイド準同型であることが自動になる。 -/
theorem cfp_frobTrivial_iff (A : CfpCat P G) :
    IsFrobeniusTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔ IsFrobeniusTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · rintro ⟨ζ, hdeg, hprop⟩
    refine ⟨⟨⟨fun n => CfpCat.fst (ζ n), ?_⟩, ?_⟩, hdeg, fun n =>
      ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ (hprop n).1,
       (cfp_frobType_iff P G hG hD' hcC hcD' _ (hprop n).2.2).mp (hprop n).2⟩⟩
    · show CfpCat.fst (ζ 1) = 𝟙 _
      rw [ζ.map_one]; rfl
    · intro m n
      show CfpCat.fst (ζ (m * n)) = CfpCat.fst (ζ n) ≫ CfpCat.fst (ζ m)
      rw [ζ.map_mul]; rfl
  · rintro ⟨ζ₀, hdeg, hprop⟩
    have hsq : ∀ n : ℕ+, P.proj.map (ζ₀ n) ≫ A.obj.hom
        = A.obj.hom ≫ G.map (𝟙 A.obj.right) := fun n => by
      rw [show P.proj.map (ζ₀ n) = 𝟙 _ from (hprop n).1.trans (P.Base_id _),
        G.map_id, Category.comp_id, Category.id_comp]
    refine ⟨⟨⟨fun n => InducedCategory.homMk ⟨ζ₀ n, 𝟙 _, hsq n⟩, ?_⟩, ?_⟩, hdeg, fun n =>
      ⟨by show CfpCat.snd _ = 𝟙 _; rfl,
       (cfp_frobType_iff P G hG hD' hcC hcD' _ (by show IsIso (𝟙 _); infer_instance)).mpr (hprop n).2⟩⟩
    · refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
      · show (ζ₀ 1 : A.obj.left ⟶ A.obj.left) = 𝟙 _
        rw [ζ₀.map_one]; rfl
      · show (𝟙 A.obj.right) = 𝟙 _; rfl
    · intro m n
      refine InducedCategory.hom_ext (CommaMorphism.ext ?_ ?_)
      · show (ζ₀ (m * n) : A.obj.left ⟶ A.obj.left) = ζ₀ n ≫ ζ₀ m
        rw [ζ₀.map_mul]; rfl
      · show (𝟙 A.obj.right) = 𝟙 _ ≫ 𝟙 _
        simp

/-- **(v)** —— **sub-quasi-Frobenius-trivial** は射影で決まる(★(a) 鎖型)。 -/
theorem cfp_subQuasiFrobTrivial_iff (A : CfpCat P G) :
    IsSubQuasiFrobeniusTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔
      IsSubQuasiFrobeniusTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · rintro ⟨Dd', α, hca, hps, hq⟩
    exact ⟨Dd'.obj.left, CfpCat.fst α, (cfp_coAngular_iff P G hG hD' hcC hcD' α).mp hca,
      (cfp_preStep_iff P G hG hD' hcC hcD' α hps.2).mp hps,
      (cfp_quasiFrobTrivial_iff P G hG hD' hcC hcD' Dd').mp hq⟩
  · rintro ⟨Dd₀, α₀, hca, hps, hq⟩
    haveI hαb : IsIso (P.proj.map α₀) := hps.2
    have hzi : IsIso (P.proj.map α₀ ≫ A.obj.hom) := inferInstance
    refine ⟨⟨⟨Dd₀, A.obj.right, P.proj.map α₀ ≫ A.obj.hom⟩, hzi⟩,
      InducedCategory.homMk ⟨α₀, 𝟙 _, by simp⟩, ?_, ⟨hps.1, by show IsIso (𝟙 _); infer_instance⟩,
      ?_⟩
    · exact cfp_coAngular_of P G hG hD' hcC hcD' _ hca
    · exact (cfp_quasiFrobTrivial_iff P G hG hD' hcC hcD' _).mpr hq

/-- **(v)** —— **unit-trivial** は射影で決まる(★(b) `𝒟'` 成分固定型)。 -/
theorem cfp_unitTrivial_iff (A : CfpCat P G) :
    IsUnitTrivial (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔ IsUnitTrivial P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  constructor
  · intro h
    refine Submonoid.eq_bot_iff_forall _ |>.mpr ?_
    intro x₀ hx₀
    have hid : P.proj.map (x₀ : A.obj.left ⟶ A.obj.left) = 𝟙 _ :=
      hx₀.1.1.trans (P.Base_id _)
    have hsq : P.proj.map (x₀ : A.obj.left ⟶ A.obj.left) ≫ A.obj.hom
        = A.obj.hom ≫ G.map (𝟙 A.obj.right) := by
      rw [hid, G.map_id, Category.comp_id, Category.id_comp]
    haveI hxi : IsIso (x₀ : A.obj.left ⟶ A.obj.left) := (isUnit_iff_isIso _).mp hx₀.2
    have hmem : (InducedCategory.homMk ⟨(x₀ : A.obj.left ⟶ A.obj.left), 𝟙 _, hsq⟩ : End A)
        ∈ OTimes (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
      refine ⟨⟨by show CfpCat.snd _ = 𝟙 _; rfl, hx₀.1.2⟩, (isUnit_iff_isIso _).mpr ?_⟩
      exact cfp_isIso_of P G _ hxi (by show IsIso (𝟙 _); infer_instance)
    have := (Submonoid.eq_bot_iff_forall _).mp h _ hmem
    exact congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) this
  · intro h
    refine Submonoid.eq_bot_iff_forall _ |>.mpr ?_
    intro x hx
    have hfst : CfpCat.fst (x : End A) ∈ OTimes P A.obj.left := by
      refine ⟨⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hx.1.1, hx.1.2⟩, (isUnit_iff_isIso _).mpr ?_⟩
      exact cfp_isIso_fst P G _ ((isUnit_iff_isIso _).mp hx.2)
    have h1 : CfpCat.fst (x : End A) = 𝟙 _ := (Submonoid.eq_bot_iff_forall _).mp h _ hfst
    exact InducedCategory.hom_ext (CommaMorphism.ext h1 hx.1.1)

/-! ### ★★**第3の基準** —— モノイドの性質として定義された型

★(a) 鎖 / (b) `𝒟'` 成分固定 は「**射が絡む定義**」用の基準である。
`group-like` のように **`Φ(A_𝒟)` というモノイドの性質**として定義された型には効かない。
そこに効くのは:

★**底の同型 `α : Base A ≅ G A'` が誘導する加法モノイド同型 `Φ(α)` に沿って移す**。

★**基準が3つに増えた理由**は「定義の形が3種類ある」ことにある:
射の存在/全称(→鎖)、自己射の条件(→`𝒟'` 成分固定)、底のモノイドの条件(→`Φ(α)`)。 -/

/-- **(v)** —— **group-like** は射影で決まる(★第3の基準)。 -/
theorem cfp_groupLike_iff (A : CfpCat P G) :
    IsGroupLikeObj (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔ IsGroupLikeObj P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  show IsGroupLike (Φ.val (G.obj A.obj.right)) ↔ IsGroupLike (Φ.val (P.proj.obj A.obj.left))
  rw [isGroupLike_iff, isGroupLike_iff]
  constructor
  · intro h a
    have := (h (Φ.map (inv A.obj.hom) a)).map (Φ.map A.obj.hom)
    rwa [← Φ.map_comp (inv A.obj.hom) A.obj.hom, IsIso.hom_inv_id, Φ.map_id] at this
  · intro h a
    have := (h (Φ.map A.obj.hom a)).map (Φ.map (inv A.obj.hom))
    rwa [← Φ.map_comp A.obj.hom (inv A.obj.hom), IsIso.inv_hom_id, Φ.map_id] at this

/-- ★`𝒞'` の自己射モノイドから `𝒞` のそれへのモノイド準同型(★`End` の積は `x * y = y ≫ x`)。 -/
def cfpEndHom (A : CfpCat P G) : End A →* End A.obj.left where
  toFun f := CfpCat.fst f
  map_one' := rfl
  map_mul' _ _ := rfl

theorem cfpEndHom_pow {A : CfpCat P G} (x : End A) (k : ℕ) :
    cfpEndHom P G A (x ^ k) = (cfpEndHom P G A x) ^ k :=
  map_pow (cfpEndHom P G A) x k

/-! ★**測定**: `CfpCat.fst (x ^ k) = (CfpCat.fst x) ^ k` を
`CfpCat.fst` の綴りで述べようとすると、型注釈を付けても
`HPow (A.obj.left ⟶ A.obj.left) ℕ` が合成できない
(`End X` の `Monoid` インスタンスは `End X` の綴りにしか付かない)。
★**したがって `cfpEndHom_pow`(返り値の型が `End A.obj.left` である形)を
`show` で当てるのが正しい**。`Frobenius-normalized` はこの形で書き直せるはずだが、未着手。 -/

/-- ★`𝒟'` 成分の方のモノイド準同型。 -/
def cfpEndHomSnd (A : CfpCat P G) : End A →* End A.obj.right where
  toFun f := CfpCat.snd f
  map_one' := rfl
  map_mul' _ _ := rfl

theorem cfpEndHomSnd_pow {A : CfpCat P G} (x : End A) (k : ℕ) :
    cfpEndHomSnd P G A (x ^ k) = (cfpEndHomSnd P G A x) ^ k :=
  map_pow (cfpEndHomSnd P G A) x k

/-! ### ★(参考) `Frobenius-normalized` の以前の記録

★**数学は自明に近い**: `𝒪^▷` の元も base-identity 自己射も `𝒟'` 成分が `𝟙` に固定されるので、
`𝒞` の等式がそのまま両成分に降りる。

★**Lean で止まる理由(特定済み)**: 原文の条件は `α ^ deg_Fr(φ)` を含み、
**`End A` と `A ⟶ A` は同じ型の2つの綴りで、`Monoid`(したがって `^`)は
`End A` の綴りにしか付かない**。`CfpCat.fst` の引数型は `X ⟶ Y` なので、
`CfpCat.fst (α ^ k)` と書くと `α ^ k` が `A ⟶ A` として推論され
`HPow (A ⟶ A) ℕ` が合成できない。

★試して駄目だったもの: 型注釈 `(… : End A)` / `cfpEndHom`(返り値型 `End A.obj.left`)経由 /
`show` で目標を書き換え / `map_mul`(`≫` は `*` ではないので不適) / `one_pow`(`𝟙` と `1` の綴り)。
★**`cfpEndHom` / `cfpEndHom_pow` / `cfpEndHomSnd` / `cfpEndHomSnd_pow` は残してある。**
次に当たる人が同じ道を繰り返さずに済む形で記録する。
-/


/-- ★`𝒞'` で base-isomorphic なら `𝒞` でも base-isomorphic(★片向き)。

★★**逆は言えない**: `𝒞` の同型 `proj A.left ≅ proj B.left` から
`A.right ≅ B.right`(`𝒟'` の同型)を作るには **`G` が同型を反映する**必要があり、
仮定にない。★これも「未解決 4 件」と同じ根である。 -/
theorem cfp_baseIsomorphic_of (A B : CfpCat P G)
    (h : BaseIsomorphic (cfpPreFrobenioid P G hG hD' hcC hcD') A B) :
    BaseIsomorphic P A.obj.left B.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  haveI hB : IsIso B.obj.hom := B.property
  obtain ⟨w⟩ := h
  exact ⟨(@asIso _ _ _ _ A.obj.hom hA) ≪≫ G.mapIso w ≪≫ (@asIso _ _ _ _ B.obj.hom hB).symm⟩

/-! ## ★(vi) —— **片向きだけ**の3タイプ

原文 (FrdI p.28):
> (vi) A object of C is Aut-ample (respectively, Autsub-ample; End-ample) if

★★**なぜ片向きなのか**が基準で説明できる:
`Aut-ample` 等は「**`𝒟` の自己射が `𝒞` から来る**」という**全射性**の条件である。
`𝒞'` について言うには `𝒟'` の自己射 `g` を `𝒟` へ送って(`G g`)`𝒞` の全射性を使えばよい ——
**`G` は関手なので送れる**。
★逆向きは「`𝒟` の自己射 `g₀` に対応する `𝒟'` の自己射」を要求するので
**`G` の充満性**が要り、仮定にない。★**原文が片向きでしか述べていないのは正しい。** -/

/-! ### ★★**#5 の原因の特定** —— 「インスタンス合成の失敗」は独立の症状ではなかった

★`IsAutAmple P' A` の `g` の型は `End ((cfpToElem …).obj A).base` であり、
**`A.obj.right` に簡約されるのは `cfpToElem` を展開したあと**である。
したがって `G.map g` の型が「もう一つの綴り」になり、
`A.obj.hom ≫ G.map g ≫ w` の `IsIso` を探すときに
文脈の `hA : IsIso A.obj.hom` と**綴りが合わず**合成が失敗する。

★★**つまり #5 は #1(2つの綴り)が「インスタンス探索」を通して現れたもの**であり、
対処も #3 と同じ ——**綴りの決まった変数を先に導入する**。
★**表は 7 行のままでよい。#5 の「原因」欄を #1 と同じに直すのが正しい。** -/

/-- ★底恒等な `𝒞` の自己射を `𝒞'` へ持ち上げる(`𝒟'` 成分は `𝟙`)。

★★**基準の6例目**: `base-identity` は `𝒟'` 成分を `𝟙` に固定するので、
持ち上げに選択の余地がなく、四角形も `Base φ₀ = 𝟙` だけで閉じる。 -/
noncomputable def cfpLiftBaseId (A : CfpCat P G) {φ₀ : End A.obj.left}
    (h : IsBaseIdentity P φ₀) : End A :=
  InducedCategory.homMk (⟨(φ₀ : A.obj.left ⟶ A.obj.left), 𝟙 A.obj.right, by
    have h1 : P.proj.map ((φ₀ : A.obj.left ⟶ A.obj.left)) = 𝟙 (P.proj.obj A.obj.left) :=
      (h : P.Base _ = P.Base (𝟙 A.obj.left)).trans (P.Base_id A.obj.left)
    show P.proj.map ((φ₀ : A.obj.left ⟶ A.obj.left)) ≫ A.obj.hom
      = A.obj.hom ≫ G.map (𝟙 A.obj.right)
    rw [h1, Category.id_comp, G.map_id, Category.comp_id]⟩ :
    (A : CfpCat P G).obj ⟶ A.obj)

@[simp] theorem cfpLiftBaseId_fst (A : CfpCat P G) {φ₀ : End A.obj.left}
    (h : IsBaseIdentity P φ₀) :
    CfpCat.fst (cfpLiftBaseId P G A h) = ((φ₀ : A.obj.left ⟶ A.obj.left)) := rfl

@[simp] theorem cfpLiftBaseId_snd (A : CfpCat P G) {φ₀ : End A.obj.left}
    (h : IsBaseIdentity P φ₀) :
    CfpCat.snd (cfpLiftBaseId P G A h) = 𝟙 A.obj.right := rfl

/-- **(v)** —— **Frobenius-normalized** は射影で決まる。

★`𝒪^▷` の元も底恒等な自己射も `𝒟'` 成分が `𝟙` なので、等式は `𝒞` 成分だけの話になる
(`𝒟'` 成分は `𝟙 ≫ 𝟙 = 𝟙 ≫ 𝟙`)。★`⟹` は `cfpLiftBaseId` で持ち上げるだけ。 -/
theorem cfp_frobNormalized_iff (A : CfpCat P G) :
    IsFrobeniusNormalized (cfpPreFrobenioid P G hG hD' hcC hcD') A
      ↔ IsFrobeniusNormalized P A.obj.left := by
  constructor
  · intro hfn φ₀ h₀ α₀ hα₀
    have hbα : IsBaseIdentity P ((α₀ : A.obj.left ⟶ A.obj.left)) := hα₀.1
    have hmem : (cfpLiftBaseId P G A hbα) ∈ OTri (cfpPreFrobenioid P G hG hD' hcC hcD') A :=
      ⟨rfl, hα₀.2⟩
    have hL := hfn (cfpLiftBaseId P G A h₀) rfl (cfpLiftBaseId P G A hbα) hmem
    have hL' : ((cfpLiftBaseId P G A hbα)
          ^ ((cfpPreFrobenioid P G hG hD' hcC hcD').degFr (cfpLiftBaseId P G A h₀) : ℕ))
        * (cfpLiftBaseId P G A h₀)
        = (cfpLiftBaseId P G A h₀) * (cfpLiftBaseId P G A hbα) := hL
    have h2 := congrArg (cfpEndHom P G A) hL'
    rw [map_mul, map_mul, map_pow] at h2
    exact h2
  · intro hfn φ hφ α hα
    have hfα : ((CfpCat.fst ((α : A ⟶ A)) : End A.obj.left)) ∈ OTri P A.obj.left :=
      ⟨cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hα.1, hα.2⟩
    have hL := hfn (CfpCat.fst ((φ : A ⟶ A)))
      (cfp_baseIdentity_fst P G hG hD' hcC hcD' _ hφ) _ hfα
    have hleft : cfpEndHom P G A
        ((α ^ ((cfpPreFrobenioid P G hG hD' hcC hcD').degFr ((φ : A ⟶ A)) : ℕ)) * φ)
        = cfpEndHom P G A (φ * α) := by
      rw [map_mul, map_mul, map_pow]
      exact hL
    have hright : cfpEndHomSnd P G A
        ((α ^ ((cfpPreFrobenioid P G hG hD' hcC hcD').degFr ((φ : A ⟶ A)) : ℕ)) * φ)
        = cfpEndHomSnd P G A (φ * α) := by
      rw [map_mul, map_mul, map_pow,
        show cfpEndHomSnd P G A φ = 1 from hφ, show cfpEndHomSnd P G A α = 1 from hα.1]
      exact congrArg (fun t => t * (1 : End A.obj.right)) (one_pow _)
    exact InducedCategory.hom_ext (CommaMorphism.ext hleft hright)

/-! ### ★(v) の `perfect` —— 対象の持ち上げ道具

★★`perfect` の主張は「`A` と底が同型な対象すべて」を量化するので、
**`𝒞` の対象を `𝒞'` へ持ち上げる**道具が要る。底が同型な射に沿えば持ち上がる。 -/

/-- ★底が同型な `𝒞` の射 `f : Z ⟶ X_𝒞` について、始域 `Z` を `𝒞'` へ持ち上げる。 -/
noncomputable def cfpLiftDom (X : CfpCat P G) {Z : C} (f : Z ⟶ X.obj.left)
    (hf : IsIso (P.proj.map f)) : CfpCat P G :=
  haveI : IsIso X.obj.hom := X.property
  haveI := hf
  have hxi : IsIso (P.proj.map f ≫ X.obj.hom) := inferInstance
  ⟨⟨Z, X.obj.right, P.proj.map f ≫ X.obj.hom⟩, hxi⟩

/-- ★その射も `𝒞'` へ持ち上がる(`𝒟'` 成分は `𝟙`)。 -/
noncomputable def cfpLiftHom (X : CfpCat P G) {Z : C} (f : Z ⟶ X.obj.left)
    (hf : IsIso (P.proj.map f)) : cfpLiftDom P G X f hf ⟶ X :=
  InducedCategory.homMk (⟨f, 𝟙 X.obj.right, by
    show P.proj.map f ≫ X.obj.hom = (P.proj.map f ≫ X.obj.hom) ≫ G.map (𝟙 X.obj.right)
    rw [G.map_id, Category.comp_id]⟩ : (cfpLiftDom P G X f hf).obj ⟶ X.obj)

@[simp] theorem cfpLiftHom_fst (X : CfpCat P G) {Z : C} (f : Z ⟶ X.obj.left)
    (hf : IsIso (P.proj.map f)) : CfpCat.fst (cfpLiftHom P G X f hf) = f := rfl

@[simp] theorem cfpLiftHom_snd (X : CfpCat P G) {Z : C} (f : Z ⟶ X.obj.left)
    (hf : IsIso (P.proj.map f)) : CfpCat.snd (cfpLiftHom P G X f hf) = 𝟙 X.obj.right := rfl

/-- ★底が同型な `𝒞` の射 `f : X_𝒞 ⟶ Z` について、終域 `Z` を `𝒞'` へ持ち上げる。 -/
noncomputable def cfpLiftCod (X : CfpCat P G) {Z : C} (f : X.obj.left ⟶ Z)
    (hf : IsIso (P.proj.map f)) : CfpCat P G :=
  haveI : IsIso X.obj.hom := X.property
  haveI := hf
  have hxi : IsIso (inv (P.proj.map f) ≫ X.obj.hom) := inferInstance
  ⟨⟨Z, X.obj.right, inv (P.proj.map f) ≫ X.obj.hom⟩, hxi⟩

/-- ★その射も `𝒞'` へ持ち上がる(`𝒟'` 成分は `𝟙`)。 -/
noncomputable def cfpLiftHomCod (X : CfpCat P G) {Z : C} (f : X.obj.left ⟶ Z)
    (hf : IsIso (P.proj.map f)) : X ⟶ cfpLiftCod P G X f hf :=
  haveI := hf
  InducedCategory.homMk (⟨f, 𝟙 X.obj.right, by
    show P.proj.map f ≫ (inv (P.proj.map f) ≫ X.obj.hom)
      = X.obj.hom ≫ G.map (𝟙 X.obj.right)
    rw [G.map_id, Category.comp_id, ← Category.assoc, IsIso.hom_inv_id,
      Category.id_comp]⟩ : X.obj ⟶ (cfpLiftCod P G X f hf).obj)

@[simp] theorem cfpLiftHomCod_fst (X : CfpCat P G) {Z : C} (f : X.obj.left ⟶ Z)
    (hf : IsIso (P.proj.map f)) : CfpCat.fst (cfpLiftHomCod P G X f hf) = f := rfl

@[simp] theorem cfpLiftHomCod_snd (X : CfpCat P G) {Z : C} (f : X.obj.left ⟶ Z)
    (hf : IsIso (P.proj.map f)) : CfpCat.snd (cfpLiftHomCod P G X f hf) = 𝟙 X.obj.right := rfl

/-- ★`𝒞` 側で決まった pre-step を `𝒞'` へ持ち上げる。

★★**`𝒟'` 成分は四角形から一意に決まる** —— `φ₂` は Frobenius 型なので
`𝒟'` 成分が同型であり、`snd ψ = snd φ₁ ≫ snd ψ' ≫ (snd φ₂)⁻¹` しかない。 -/
noncomputable def cfpPerfectLift {B₁ B₁' B₂ B₂' : CfpCat P G} (φ₁ : B₁ ⟶ B₁') (φ₂ : B₂ ⟶ B₂')
    [hs2 : IsIso (CfpCat.snd φ₂)]
    (ψ' : B₁' ⟶ B₂') (ψ₀ : B₁.obj.left ⟶ B₂.obj.left)
    (hψ₀e : CfpCat.fst φ₁ ≫ CfpCat.fst ψ' = ψ₀ ≫ CfpCat.fst φ₂) : B₁ ⟶ B₂ :=
  InducedCategory.homMk (⟨ψ₀, CfpCat.snd φ₁ ≫ CfpCat.snd ψ' ≫ inv (CfpCat.snd φ₂), by
    haveI hg2 : IsIso (G.map (CfpCat.snd φ₂)) := inferInstance
    refine (cancel_mono (G.map (CfpCat.snd φ₂))).mp ?_
    have hL : (P.proj.map ψ₀ ≫ B₂.obj.hom) ≫ G.map (CfpCat.snd φ₂)
        = P.proj.map (CfpCat.fst φ₁) ≫ P.proj.map (CfpCat.fst ψ') ≫ B₂'.obj.hom := by
      rw [Category.assoc, ← cfp_square φ₂, ← Category.assoc, ← P.proj.map_comp, ← hψ₀e,
        P.proj.map_comp, Category.assoc]
    rw [hL, cfp_square ψ', ← Category.assoc, cfp_square φ₁, Category.assoc, ← G.map_comp,
      Category.assoc, ← G.map_comp, Category.assoc, Category.assoc, IsIso.inv_hom_id,
      Category.comp_id]⟩ : B₁.obj ⟶ B₂.obj)

@[simp] theorem cfpPerfectLift_fst {B₁ B₁' B₂ B₂' : CfpCat P G} (φ₁ : B₁ ⟶ B₁') (φ₂ : B₂ ⟶ B₂')
    [IsIso (CfpCat.snd φ₂)] (ψ' : B₁' ⟶ B₂') (ψ₀ : B₁.obj.left ⟶ B₂.obj.left)
    (hψ₀e : CfpCat.fst φ₁ ≫ CfpCat.fst ψ' = ψ₀ ≫ CfpCat.fst φ₂) :
    CfpCat.fst (cfpPerfectLift P G φ₁ φ₂ ψ' ψ₀ hψ₀e) = ψ₀ := rfl

@[simp] theorem cfpPerfectLift_snd {B₁ B₁' B₂ B₂' : CfpCat P G} (φ₁ : B₁ ⟶ B₁') (φ₂ : B₂ ⟶ B₂')
    [IsIso (CfpCat.snd φ₂)] (ψ' : B₁' ⟶ B₂') (ψ₀ : B₁.obj.left ⟶ B₂.obj.left)
    (hψ₀e : CfpCat.fst φ₁ ≫ CfpCat.fst ψ' = ψ₀ ≫ CfpCat.fst φ₂) :
    CfpCat.snd (cfpPerfectLift P G φ₁ φ₂ ψ' ψ₀ hψ₀e)
      = CfpCat.snd φ₁ ≫ CfpCat.snd ψ' ≫ inv (CfpCat.snd φ₂) := rfl

/-- **(v)** —— **perfect** は射影から降りる。 -/
theorem cfp_perfect_mpr (A : CfpCat P G) (h : IsPerfectObj P A.obj.left) :
    IsPerfectObj (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
  intro n
  refine ⟨?_, ?_⟩
  · intro B hB
    obtain ⟨B₀, φ, hφ, hdeg⟩ := (h n).1 B.obj.left
      (cfp_baseIsomorphic_of P G hG hD' hcC hcD' A B hB)
    refine ⟨cfpLiftDom P G B φ hφ.2, cfpLiftHom P G B φ hφ.2, ?_, hdeg⟩
    exact (cfp_frobType_iff P G hG hD' hcC hcD' _
      (by show IsIso (𝟙 B.obj.right); infer_instance)).mpr hφ
  · intro B₁ B₁' B₂ B₂' φ₁ φ₂ hφ₁ hd₁ hφ₂ hd₂ hb₁ hb₂ ψ' hψ'
    haveI hs2 : IsIso (CfpCat.snd φ₂) := hφ₂.2
    haveI hs1 : IsIso (CfpCat.snd φ₁) := hφ₁.2
    haveI hsp : IsIso (CfpCat.snd ψ') := hψ'.2
    obtain ⟨ψ₀, ⟨hψ₀s, hψ₀e⟩, hψ₀u⟩ := (h n).2 B₁.obj.left B₁'.obj.left B₂.obj.left B₂'.obj.left
      (CfpCat.fst φ₁) (CfpCat.fst φ₂)
      ((cfp_frobType_iff P G hG hD' hcC hcD' φ₁ hφ₁.2).mp hφ₁) hd₁
      ((cfp_frobType_iff P G hG hD' hcC hcD' φ₂ hφ₂.2).mp hφ₂) hd₂
      (cfp_baseIsomorphic_of P G hG hD' hcC hcD' B₁ A hb₁)
      (cfp_baseIsomorphic_of P G hG hD' hcC hcD' B₂ A hb₂)
      (CfpCat.fst ψ') ((cfp_preStep_iff P G hG hD' hcC hcD' ψ' hψ'.2).mp hψ')
    refine ⟨cfpPerfectLift P G φ₁ φ₂ ψ' ψ₀ hψ₀e, ⟨⟨hψ₀s.1, ?_⟩, ?_⟩, ?_⟩
    · show IsIso (CfpCat.snd φ₁ ≫ CfpCat.snd ψ' ≫ inv (CfpCat.snd φ₂))
      infer_instance
    · refine InducedCategory.hom_ext (CommaMorphism.ext hψ₀e ?_)
      show CfpCat.snd φ₁ ≫ CfpCat.snd ψ'
        = (CfpCat.snd φ₁ ≫ CfpCat.snd ψ' ≫ inv (CfpCat.snd φ₂)) ≫ CfpCat.snd φ₂
      rw [Category.assoc, Category.assoc, IsIso.inv_hom_id, Category.comp_id]
    · intro ψ₂ hψ₂
      have hf : CfpCat.fst ψ₂ = ψ₀ :=
        hψ₀u (CfpCat.fst ψ₂)
          ⟨(cfp_preStep_iff P G hG hD' hcC hcD' ψ₂ hψ₂.1.2).mp hψ₂.1,
            congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hψ₂.2⟩
      have hs : CfpCat.snd ψ₂ = CfpCat.snd φ₁ ≫ CfpCat.snd ψ' ≫ inv (CfpCat.snd φ₂) := by
        have h2 : CfpCat.snd φ₁ ≫ CfpCat.snd ψ' = CfpCat.snd ψ₂ ≫ CfpCat.snd φ₂ :=
          congrArg (fun t => CommaMorphism.right (InducedCategory.Hom.hom t)) hψ₂.2
        rw [← Category.assoc, h2, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
      exact InducedCategory.hom_ext (CommaMorphism.ext hf hs)

/-- **(v)** —— **perfect** は射影へ上がる。

★★`𝒞` の対象はすべて `𝒞'` へ持ち上がる —— 底の同型を `𝒟'` 成分 `A_𝒟'` へ運ぶだけ。
★連鎖 `B₁ →(φ₁) B₁' →(ψ') B₂' ←(φ₂) B₂` は**すべて底同型**なので、
`B₁'` を 1 つ持ち上げれば残りは `cfpLiftDom` / `cfpLiftCod` で決まり、
**`𝒟'` 成分はすべて `𝟙`** になる。 -/
theorem cfp_perfect_mp (A : CfpCat P G)
    (h : IsPerfectObj (cfpPreFrobenioid P G hG hD' hcC hcD') A) :
    IsPerfectObj P A.obj.left := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨a, ha⟩ : ∃ a : (P.toElem.obj A.obj.left).base ⟶ G.obj A.obj.right,
      a = A.obj.hom := ⟨A.obj.hom, rfl⟩
  haveI hai : IsIso a := ha ▸ hA
  intro n
  refine ⟨?_, ?_⟩
  · intro B hB
    obtain ⟨e⟩ := hB
    haveI hwi : IsIso (e.inv ≫ a) := inferInstance
    obtain ⟨B₀', φ', hφ', hdeg⟩ := (h n).1 (cfpMk P G B A.obj.right (e.inv ≫ a)) ⟨Iso.refl _⟩
    exact ⟨B₀'.obj.left, CfpCat.fst φ',
      (cfp_frobType_iff P G hG hD' hcC hcD' φ' hφ'.2).mp hφ', hdeg⟩
  · intro B₁ B₁' B₂ B₂' φ₁ φ₂ hφ₁ hd₁ hφ₂ hd₂ hb₁ hb₂ ψ' hψ'
    haveI hp1 : IsIso (P.proj.map φ₁) := hφ₁.2
    haveI hp2 : IsIso (P.proj.map φ₂) := hφ₂.2
    haveI hpp : IsIso (P.proj.map ψ') := hψ'.2
    obtain ⟨e₁⟩ := hb₁
    obtain ⟨b, hb⟩ : ∃ b : P.proj.obj B₁ ⟶ G.obj A.obj.right, b = e₁.hom ≫ a := ⟨_, rfl⟩
    haveI hbi : IsIso b := by
      rw [hb]; haveI : IsIso (e₁.hom ≫ a) := inferInstance; exact this
    haveI hw1 : IsIso (inv (P.proj.map φ₁) ≫ b) := inferInstance
    set L1' : CfpCat P G := cfpMk P G B₁' A.obj.right (inv (P.proj.map φ₁) ≫ b) with hL1'
    set L2' : CfpCat P G := cfpLiftCod P G L1' ψ' hpp with hL2'
    set L1 : CfpCat P G := cfpLiftDom P G L1' φ₁ hp1 with hL1
    set L2 : CfpCat P G := cfpLiftDom P G L2' φ₂ hp2 with hL2
    obtain ⟨Ψ, ⟨hΨs, hΨe⟩, hΨu⟩ := (h n).2 L1 L1' L2 L2'
      (cfpLiftHom P G L1' φ₁ hp1) (cfpLiftHom P G L2' φ₂ hp2)
      ((cfp_frobType_iff P G hG hD' hcC hcD' _
        (by show IsIso (𝟙 A.obj.right); infer_instance)).mpr hφ₁) hd₁
      ((cfp_frobType_iff P G hG hD' hcC hcD' _
        (by show IsIso (𝟙 A.obj.right); infer_instance)).mpr hφ₂) hd₂
      ⟨Iso.refl _⟩ ⟨Iso.refl _⟩
      (cfpLiftHomCod P G L1' ψ' hpp)
      ⟨hψ'.1, by show IsIso (𝟙 A.obj.right); infer_instance⟩
    refine ⟨CfpCat.fst Ψ, ⟨(cfp_preStep_iff P G hG hD' hcC hcD' Ψ hΨs.2).mp hΨs, ?_⟩, ?_⟩
    · exact congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t)) hΨe
    · intro ψ₂ hψ₂
      have hcomm : P.proj.map (ψ₂ ≫ φ₂) = P.proj.map (φ₁ ≫ ψ') :=
        congrArg P.proj.map hψ₂.2.symm
      have hsq : P.proj.map ψ₂ ≫ L2.obj.hom = L1.obj.hom ≫ G.map (𝟙 A.obj.right) := by
        show P.proj.map ψ₂ ≫ (P.proj.map φ₂ ≫ (inv (P.proj.map ψ')
            ≫ (inv (P.proj.map φ₁) ≫ b)))
          = (P.proj.map φ₁ ≫ (inv (P.proj.map φ₁) ≫ b)) ≫ G.map (𝟙 A.obj.right)
        rw [G.map_id, Category.comp_id, ← Category.assoc, ← P.proj.map_comp, hcomm,
          P.proj.map_comp]
        simp only [Category.assoc]
        rw [← Category.assoc (P.proj.map ψ') (inv (P.proj.map ψ')), IsIso.hom_inv_id,
          Category.id_comp]
      refine congrArg (fun t => CommaMorphism.left (InducedCategory.Hom.hom t))
        (hΨu (cfpHom P G ψ₂ (𝟙 A.obj.right) hsq) ⟨⟨hψ₂.1.1, ?_⟩, ?_⟩)
      · show IsIso (𝟙 A.obj.right)
        infer_instance
      · exact InducedCategory.hom_ext (CommaMorphism.ext hψ₂.2 rfl)

/-- **(v)** —— **perfect** は射影で決まる(両向き)。 -/
theorem cfp_perfect_iff (A : CfpCat P G) :
    IsPerfectObj (cfpPreFrobenioid P G hG hD' hcC hcD') A ↔ IsPerfectObj P A.obj.left :=
  ⟨cfp_perfect_mp P G hG hD' hcC hcD' A, cfp_perfect_mpr P G hG hD' hcC hcD' A⟩

/-- **(vi)** —— **End-ample** は射影から降りる(★片向き)。 -/
theorem cfp_endAmple_of (A : CfpCat P G) (h : IsEndAmple P A.obj.left) :
    IsEndAmple (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨w, hw1, hw2⟩ := hA.out
  intro g0
  -- ★#3: 綴りの決まった変数を先に導入する
  obtain ⟨g, rfl⟩ : ∃ g : End A.obj.right, g = g0 := ⟨g0, rfl⟩
  obtain ⟨φ₀, hφ₀⟩ := h (A.obj.hom ≫ G.map g ≫ w)
  refine ⟨InducedCategory.homMk ⟨φ₀, g, ?_⟩, rfl⟩
  show P.proj.map φ₀ ≫ A.obj.hom = A.obj.hom ≫ G.map g
  rw [show P.proj.map φ₀ = A.obj.hom ≫ G.map g ≫ w from hφ₀, Category.assoc,
    Category.assoc, hw2, Category.comp_id]

/-- **(vi)** —— **Aut-ample** は射影から降りる(★片向き)。 -/
theorem cfp_autAmple_of (A : CfpCat P G) (h : IsAutAmple P A.obj.left) :
    IsAutAmple (cfpPreFrobenioid P G hG hD' hcC hcD') A := by
  haveI hA : IsIso A.obj.hom := A.property
  obtain ⟨w, hw1, hw2⟩ := hA.out
  haveI hwi : IsIso w := ⟨A.obj.hom, hw2, hw1⟩
  intro g0 hg0
  obtain ⟨g, rfl⟩ : ∃ g : End A.obj.right, g = g0 := ⟨g0, rfl⟩
  haveI hgi : IsIso g := hg0
  haveI hGg : IsIso (G.map g) := inferInstance
  haveI hcomp : IsIso (A.obj.hom ≫ G.map g ≫ w) := inferInstance
  obtain ⟨φ₀, hiso, hφ₀⟩ := h (A.obj.hom ≫ G.map g ≫ w) hcomp
  refine ⟨InducedCategory.homMk ⟨φ₀, g, ?_⟩, ?_, rfl⟩
  · show P.proj.map φ₀ ≫ A.obj.hom = A.obj.hom ≫ G.map g
    rw [show P.proj.map φ₀ = A.obj.hom ≫ G.map g ≫ w from hφ₀, Category.assoc,
      Category.assoc, hw2, Category.comp_id]
  · exact cfp_isIso_of P G _ hiso hg0

/-! ## ★(ii) —— `Definition 1.3` の 21 条の移送

原文 (FrdI p.28):
> equivalences, the conditions of Definition 1.3 follow via a routine verification. Thus,

★**机上の仕分けを実装で検証する段**である。まず辞書から直に出る条から。 -/

end Dict
end ABC3.Found.FrdI
