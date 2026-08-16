import ABC3.Found.FrdI.Prop15

/-!
# 捻れ積 `𝔽_Φ ⋉ G` —— `[FrdI] Proposition 1.14, (iii)` の反例の候補

★`Gap/FrdI/Section1.lean` の `Gap_1_14_iii` を ③(`sourceGap`)に上げるための構成。

## ★なぜ単なる直積では駄目か

合成を `(f, g) ≫ (f′, g′) = (f ≫ f′, g · g′)` と定めると、
`Frobenius-normalized` の等式 `ζ ≫ α^d = α ≫ ζ` の `G` 成分は
`g_ζ · g_α^d = g_α · g_ζ` となり、★**`g_α^{d−1} = 1` をすべての `d` について
強制する**ので `G` が自明になってしまう。

## ★★捻れ積

合成を **Frobenius 次数で捻る**:
```
(f, g) ≫ (f′, g′) = (f ≫ f′, g ^ degFr(f′) · g′)
```
`Base`・`Div`・`degFr` は `G` 成分を無視する(射影関手 `TwProj`)。

★**`Frobenius-normalized` が自動で成り立つ**:
`α = (𝟙, a)`(次数 1)、`ζ = (z, 1)`(次数 `d`)に対し
`ζ ≫ α^d = (z, 1^1 · a^d) = (z, a^d)`、`α ≫ ζ = (z, a^d · 1) = (z, a^d)`。

★★**そして `ζ` は mono でない**: `(f₀, u) ≫ (z, 1) = (f₀ ≫ z, u^d)` なので、
`u^d = v^d` から `u = v` は出ない。★**`G` が `d`-捻れを持てば破れる。**

## ★残る仕事

★**この捻れ積が `Definition 1.3` の全条件を満たすことの検証**である。
`𝔽_Φ` が Frobenioid であることは `Proposition 1.5`(`elemFrob_isFrobenioid`)で
証明済みなので、各条件に `G` 成分の処理を足す形になると見ている(★未検証)。
-/

namespace ABC3.Check.FrdI

open CategoryTheory ABC3.Found.FrdI

universe v u w t

variable {D : Type u} [Category.{v} D] {Φ : MonoidOn.{v, u, w} D}
  {G : Type t} [CommGroup G]

/-- ★**捻れ積の対象** —— `𝔽_Φ` の対象そのもの(型シノニムとして包む)。

★`G` は本体に現れないが、**圏構造が `G` に依存する**ので型に含める。 -/
structure TwObj (Φ : MonoidOn.{v, u, w} D) (G : Type t) [CommGroup G] where
  /-- 元の `𝔽_Φ` の対象 -/
  ofElem : ElemFrobCat Φ

/-- ★**捻れ積の射** —— `𝔽_Φ` の射に `G` の元を添えたもの。 -/
@[ext]
structure TwHom (A B : TwObj Φ G) where
  /-- `𝔽_Φ` 成分 -/
  hom : A.ofElem ⟶ B.ofElem
  /-- `G` 成分 -/
  unit : G

/-- ★★**捻れ積の圏構造** —— 合成の `G` 成分を **Frobenius 次数で捻る**。

★結合律は `(g ≫ h).deg = h.deg * g.deg`(`𝔽_Φ` では `rfl`)と
`G` の可換性から出る。 -/
instance twCategory : Category (TwObj Φ G) where
  Hom A B := TwHom A B
  id A := ⟨𝟙 A.ofElem, 1⟩
  comp f g := ⟨f.hom ≫ g.hom, f.unit ^ ((g.hom.deg : ℕ+) : ℕ) * g.unit⟩
  id_comp f := by
    refine TwHom.ext (Category.id_comp _) ?_
    show (1 : G) ^ ((f.hom.deg : ℕ+) : ℕ) * f.unit = f.unit
    rw [one_pow, one_mul]
  comp_id := fun {_ Y} f => by
    refine TwHom.ext (Category.comp_id _) ?_
    show f.unit ^ ((ElemFrobCat.Hom.deg (𝟙 Y.ofElem) : ℕ+) : ℕ) * (1 : G) = f.unit
    rw [show ElemFrobCat.Hom.deg (𝟙 Y.ofElem) = 1 from rfl]
    simp
  assoc f g h := by
    refine TwHom.ext (Category.assoc _ _ _) ?_
    show (f.unit ^ ((g.hom.deg : ℕ+) : ℕ) * g.unit) ^ ((h.hom.deg : ℕ+) : ℕ) * h.unit
      = f.unit ^ (((g.hom ≫ h.hom).deg : ℕ+) : ℕ)
          * (g.unit ^ ((h.hom.deg : ℕ+) : ℕ) * h.unit)
    rw [show (g.hom ≫ h.hom).deg = h.hom.deg * g.hom.deg from rfl, PNat.mul_coe,
      mul_pow, ← pow_mul, mul_comm ((h.hom.deg : ℕ+) : ℕ), mul_assoc]

@[simp] theorem twComp_hom {A B E : TwObj Φ G} (f : A ⟶ B) (g : B ⟶ E) :
    (f ≫ g).hom = f.hom ≫ g.hom := rfl

@[simp] theorem twComp_unit {A B E : TwObj Φ G} (f : A ⟶ B) (g : B ⟶ E) :
    (f ≫ g).unit = f.unit ^ ((g.hom.deg : ℕ+) : ℕ) * g.unit := rfl

@[simp] theorem twId_hom (A : TwObj Φ G) : (𝟙 A : A ⟶ A).hom = 𝟙 A.ofElem := rfl

@[simp] theorem twId_unit (A : TwObj Φ G) : (𝟙 A : A ⟶ A).unit = 1 := rfl

/-- ★**射影関手 `𝔽_Φ ⋉ G → 𝔽_Φ`** —— `G` 成分を落とす。

★★**これが pre-Frobenioid 構造を与える** —— `Base`・`Div`・`degFr` は
すべて `G` 成分を無視する。 -/
def twProj (Φ : MonoidOn.{v, u, w} D) (G : Type t) [CommGroup G] :
    TwObj Φ G ⥤ ElemFrobCat Φ where
  obj A := A.ofElem
  map f := f.hom

/-! ## ★★核心 —— 次数 `d` の射は `G` の `d`-捻れを潰す -/

/-- ★★**捻れ積では、次数 `d` の射での右合成が `G` 成分を `d` 乗する**。 -/
theorem twComp_unit_pow {A B E : TwObj Φ G} (f : A ⟶ B) (ζ : B ⟶ E)
    (hζ : ζ.unit = 1) : (f ≫ ζ).unit = f.unit ^ ((ζ.hom.deg : ℕ+) : ℕ) := by
  rw [twComp_unit, hζ, mul_one]

/-- ★★★**次数 `d` の射は mono でない** —— `G` に `d`-捻れがあるとき。

★`u ≠ 1` かつ `u ^ d = 1` なる `u ∈ G` があれば、
`(𝟙, 1) ≫ ζ = (𝟙, u) ≫ ζ` だが `(𝟙, 1) ≠ (𝟙, u)`。

★★**これが `Proposition 1.14, (iii)` の `⟸` を破る機構である** ——
そこは「prime-Frobenius 射が FSMI(したがって mono)」を要求する。 -/
theorem not_mono_twist {B E : TwObj Φ G} (ζ : B ⟶ E) (hζ : ζ.unit = 1)
    (u : G) (hu : u ≠ 1) (hup : u ^ ((ζ.hom.deg : ℕ+) : ℕ) = 1) : ¬ Mono ζ := by
  intro hm
  have hsq : (⟨𝟙 B.ofElem, 1⟩ : TwHom B B) ≫ ζ = (⟨𝟙 B.ofElem, u⟩ : TwHom B B) ≫ ζ := by
    refine TwHom.ext rfl ?_
    rw [twComp_unit, twComp_unit, hζ, mul_one, mul_one, one_pow, hup]
  have := (cancel_mono ζ).mp hsq
  exact hu (congrArg TwHom.unit this).symm

/-! ## ★捻れ積の `PreFrobenioid` 構造

★`𝔽_Φ` の側の結果に **`G` 成分の処理を足す**だけで出る。
-/

/-- ★**捻れ積も totally epimorphic**。

★`𝔽_Φ` 成分は `isTotallyEpimorphic_elemFrobCat`、
★**`G` 成分は群の消約**(指数が一致するので左から消せる)。 -/
theorem twTotallyEpimorphic (hD : IsTotallyEpimorphic D)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A)) :
    IsTotallyEpimorphic (TwObj Φ G) := by
  intro A B f
  haveI : Epi f.hom := isTotallyEpimorphic_elemFrobCat hD hint _ _ f.hom
  refine ⟨fun {Z} h₁ h₂ hh => ?_⟩
  have hhom : h₁.hom = h₂.hom := (cancel_epi f.hom).mp (congrArg TwHom.hom hh)
  refine TwHom.ext hhom ?_
  have hu := congrArg TwHom.unit hh
  rw [twComp_unit, twComp_unit, hhom] at hu
  exact mul_left_cancel hu

/-- ★**捻れ積も connected** —— `𝔽_Φ` の証明をそのまま写す。 -/
theorem twIsConnected [IsConnected D] : IsConnected (TwObj Φ G) := by
  obtain ⟨d₀⟩ := (inferInstance : Nonempty D)
  refine IsConnected.of_induct (j₀ := (⟨⟨d₀⟩⟩ : TwObj Φ G)) ?_
  intro p hp0 hstep A
  have key : ∀ d : D, (⟨⟨d⟩⟩ : TwObj Φ G) ∈ p :=
    induct_on_objects (J := D) {d | (⟨⟨d⟩⟩ : TwObj Φ G) ∈ p} hp0
      (fun {d₁ d₂} f => hstep
        (⟨(⟨f, 0, 1⟩ : (⟨d₁⟩ : ElemFrobCat Φ) ⟶ (⟨d₂⟩ : ElemFrobCat Φ)), 1⟩ :
          (⟨⟨d₁⟩⟩ : TwObj Φ G) ⟶ ⟨⟨d₂⟩⟩))
  exact key A.ofElem.base

/-- ★★**捻れ積の pre-Frobenioid 構造**。

★`𝔽_Φ` の pre-Frobenioid 構造に射影関手を前合成するだけ ——
★**`Base`・`Div`・`degFr` はすべて `G` 成分を無視する。** -/
def twPreFrobenioid (hD : IsTotallyEpimorphic D) [IsConnected D]
    (hpd : ∀ A : D, IsPreDivisorial (Φ.val A)) :
    PreFrobenioid (TwObj Φ G) Φ.charOn where
  toElem := twProj Φ G ⋙ elemFrobToChar Φ
  divisorial := Φ.charOn_isDivisorialOn hpd
  totEpiC := twTotallyEpimorphic hD (fun A => (hpd A).1)
  totEpiD := hD
  connectedC := twIsConnected
  connectedD := inferInstance

/-! ## ★同型は `𝔽_Φ` 成分だけで決まる

★★**これが「捻れ積も isotropic 型である」ことの土台**である ——
`Div = 0` の pre-step は `𝔽_Φ` で同型になり、`G` 成分は自動で逆が取れる。
-/

/-- ★★**`𝔽_Φ` 成分が同型なら、捻れ積でも同型**。

★逆射の `G` 成分は `f.unit⁻¹`。★**同型の次数は 1**(`ElemFrobCat.isIso_iff`)
なので指数が消え、両側の等式がそろう。 -/
theorem twIsIso_of {A B : TwObj Φ G} (f : A ⟶ B) (h : IsIso f.hom) : IsIso f := by
  haveI := h
  obtain ⟨-, -, hdeg⟩ := (ElemFrobCat.isIso_iff f.hom).mp h
  have hdegi : ElemFrobCat.Hom.deg (inv f.hom) = 1 := by
    have h1 := congrArg ElemFrobCat.Hom.deg (IsIso.inv_hom_id f.hom)
    rw [ElemFrobCat.degFr_comp, hdeg, one_mul] at h1
    exact h1
  refine ⟨⟨⟨inv f.hom, f.unit⁻¹⟩, ?_, ?_⟩⟩
  · refine TwHom.ext (IsIso.hom_inv_id f.hom) ?_
    show f.unit ^ ((ElemFrobCat.Hom.deg (inv f.hom) : ℕ+) : ℕ) * f.unit⁻¹ = (1 : G)
    rw [hdegi]
    simp
  · refine TwHom.ext (IsIso.inv_hom_id f.hom) ?_
    show f.unit⁻¹ ^ ((ElemFrobCat.Hom.deg f.hom : ℕ+) : ℕ) * f.unit = (1 : G)
    rw [hdeg]
    simp

/-! ## ★不変量は `𝔽_Φ` 成分だけで決まる(橋渡し)

★★**これが「`G` 成分を無視する」の正確な意味**である。 -/

section Bridge

variable (hD : IsTotallyEpimorphic D) [IsConnected D]
  (hpd : ∀ A : D, IsPreDivisorial (Φ.val A))

@[simp] theorem twBase {A B : TwObj Φ G} (f : A ⟶ B) :
    (twPreFrobenioid (G := G) hD hpd).Base f = f.hom.base := rfl

@[simp] theorem twDegFr {A B : TwObj Φ G} (f : A ⟶ B) :
    (twPreFrobenioid (G := G) hD hpd).degFr f = f.hom.deg := rfl

@[simp] theorem twDiv {A B : TwObj Φ G} (f : A ⟶ B) :
    (twPreFrobenioid (G := G) hD hpd).Div f = toChar f.hom.div := rfl

/-- ★★**捻れ積も isotropic 型** —— `Div = 0` の pre-step は同型。

★`𝔽_Φ` の側(`elemFrob_isotropic`)で `hom` が同型になり、
`twIsIso_of` が `G` 成分の逆を自動で作る。 -/
theorem twIsotropic (A : TwObj Φ G) : IsIsotropic (twPreFrobenioid (G := G) hD hpd) A := by
  intro Dd φ hisom hstep
  refine twIsIso_of φ ?_
  exact elemFrob_isotropic Φ hD hpd A.ofElem Dd.ofElem φ.hom hisom hstep

/-- ★★**捻れ積では全射が co-angular**(`Proposition 1.4, (i)`)。

★★**これで `Definition 1.3` の (iii)(a)(b) が自明になる。** -/
theorem twCoAngular {A B : TwObj Φ G} (φ : A ⟶ B) :
    IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ :=
  prop_1_4_i _ φ (fun Y _ => twIsotropic hD hpd Y)

/-- ★**捻れ積でも pre-step は mono**(`Definition 1.3, (v), (a)`)。

★`𝔽_Φ` 成分は `elemFrob_preStepMono`、★**`G` 成分は次数 1 なので単純に消約**。 -/
theorem twPreStepMono {A B : TwObj Φ G} (φ : A ⟶ B)
    (hφ : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ) : Mono φ := by
  haveI : Mono φ.hom := elemFrob_preStepMono Φ hD hpd φ.hom hφ
  refine ⟨fun {Z} f g hh => ?_⟩
  have hhom : f.hom = g.hom := (cancel_mono φ.hom).mp (congrArg TwHom.hom hh)
  refine TwHom.ext hhom ?_
  have hu := congrArg TwHom.unit hh
  rw [twComp_unit, twComp_unit, show ElemFrobCat.Hom.deg φ.hom = 1 from hφ.1] at hu
  simpa using mul_right_cancel hu

/-- ★**各次数の Frobenius 型自己射がある**(`Definition 1.3, (ii)` の存在部分)。

★`(⟨𝟙, 0, n⟩, 1)` を取ればよい。★**co-angular は `twCoAngular` が無条件に与える。** -/
theorem twFrobDegSurj (A : TwObj Φ G) (n : ℕ+) :
    ∃ (B : TwObj Φ G) (φ : A ⟶ B),
      IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) φ ∧
      (twPreFrobenioid (G := G) hD hpd).degFr φ = n := by
  refine ⟨A, ⟨⟨𝟙 A.ofElem.base, 0, n⟩, 1⟩, ⟨⟨twCoAngular hD hpd _, ?_⟩, ?_⟩, rfl⟩
  · show toChar (0 : Φ.val A.ofElem.base) = 0
    simp
  · show IsIso (𝟙 A.ofElem.base)
    infer_instance

/-- ★**(iii)(a)** co-angular は合成で閉じる —— ★**全射が co-angular なので自明**。 -/
theorem twCoAngularComp {A B E : TwObj Φ G} (ψ : A ⟶ B) (φ : B ⟶ E)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) ψ)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ) :
    IsCoAngular (twPreFrobenioid (G := G) hD hpd) (ψ ≫ φ) :=
  twCoAngular hD hpd _

/-- ★**(vii)(b)** isotropic は下流へ伝播 —— ★**全対象が isotropic なので自明**。 -/
theorem twIsotropicClosed {A B : TwObj Φ G} (_ : A ⟶ B)
    (_ : IsIsotropic (twPreFrobenioid (G := G) hD hpd) A) :
    IsIsotropic (twPreFrobenioid (G := G) hD hpd) B :=
  twIsotropic hD hpd B

/-- ★**(ii) の本質的一意性** —— 同じ次数の Frobenius 型射は同型を除いて一意。

★`𝔽_Φ` 側の `β₀` に、★**`G` 成分を `φ.unit⁻¹ * ψ.unit` と取れば四角形が閉じる**
(`β₀` は同型なので次数 1)。 -/
theorem twFrobType_hom {A B : TwObj Φ G} (φ : A ⟶ B)
    (hφ : IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) φ) :
    IsFrobeniusType (elemPreFrobenioid Φ hD hpd) φ.hom :=
  ⟨⟨prop_1_4_i _ φ.hom (fun Y _ => elemFrob_isotropic Φ hD hpd Y), hφ.1.2⟩, hφ.2⟩

theorem twFrobDegUniq (A B E : TwObj Φ G) (φ : A ⟶ B) (ψ : A ⟶ E)
    (hφ : IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) φ)
    (hψ : IsFrobeniusType (twPreFrobenioid (G := G) hD hpd) ψ)
    (hd : (twPreFrobenioid (G := G) hD hpd).degFr φ
        = (twPreFrobenioid (G := G) hD hpd).degFr ψ) :
    ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ := by
  obtain ⟨β₀, hβ₀iso, hβ₀⟩ :=
    elemFrob_frobDegUniq Φ hD hpd A.ofElem B.ofElem E.ofElem φ.hom ψ.hom
      (twFrobType_hom hD hpd φ hφ) (twFrobType_hom hD hpd ψ hψ) hd
  haveI := hβ₀iso
  obtain ⟨-, -, hdβ₀⟩ := (ElemFrobCat.isIso_iff β₀).mp hβ₀iso
  refine ⟨⟨β₀, φ.unit⁻¹ * ψ.unit⟩, twIsIso_of _ hβ₀iso, ?_⟩
  refine TwHom.ext hβ₀ ?_
  show φ.unit ^ ((ElemFrobCat.Hom.deg β₀ : ℕ+) : ℕ) * (φ.unit⁻¹ * ψ.unit) = ψ.unit
  rw [hdβ₀]
  simp

/-- ★`ζ n := (⟨𝟙, 0, n⟩, 1)` —— Frobenius-trivial 性を与えるモノイド準同型。

★★**`End` の乗法は `x * y = y ≫ x`** なので、`map_mul'` は
`ζ (n*m) = ζ m ≫ ζ n` を示すことになる。 -/
def twZeta (Y : D) : ℕ+ →* End (⟨⟨Y⟩⟩ : TwObj Φ G) where
  toFun n := ⟨⟨𝟙 Y, 0, n⟩, 1⟩
  map_one' := TwHom.ext (ElemFrobCat.Hom.ext rfl rfl rfl) rfl
  map_mul' n m := by
    refine TwHom.ext ?_ ?_
    · refine ElemFrobCat.Hom.ext ?_ ?_ rfl
      · show 𝟙 Y = 𝟙 Y ≫ 𝟙 Y
        simp
      · show (0 : Φ.val Y) = Φ.map (𝟙 Y) 0 + ((n : ℕ+) : ℕ) • (0 : Φ.val Y)
        simp
    · show (1 : G) = (1 : G) ^ ((n : ℕ+) : ℕ) * 1
      simp

/-- ★**(i)(a)** 各底に Frobenius-trivial 対象がある。 -/
theorem twBaseSurj (Y : D) :
    ∃ A : TwObj Φ G, IsFrobeniusTrivial (twPreFrobenioid (G := G) hD hpd) A ∧
      Nonempty (((twPreFrobenioid (G := G) hD hpd).toElem.obj A).base ≅ Y) := by
  refine ⟨⟨⟨Y⟩⟩, ⟨twZeta (G := G) Y, fun n => rfl, fun n => ?_⟩, ⟨Iso.refl _⟩⟩
  exact ⟨rfl, ⟨twCoAngular hD hpd _, by show toChar (0 : Φ.val Y) = 0; simp⟩,
    by show IsIso (𝟙 Y); infer_instance⟩

/-- ★**(iii)(b)** —— ★**全射が co-angular なので自明**。 -/
theorem twCoAngularOfPreStep {A' A : TwObj Φ G} (α : A' ⟶ A)
    (_ : IsCoAngular (twPreFrobenioid (G := G) hD hpd) α)
    (_ : IsPreStep (twPreFrobenioid (G := G) hD hpd) α) (φ : A' ⟶ A) :
    IsCoAngular (twPreFrobenioid (G := G) hD hpd) φ :=
  twCoAngular hD hpd φ

/-- ★**(vii)(a)** isotropic hull は `𝟙` —— ★**全対象が isotropic だから**。 -/
theorem twIsotropicHullExists (A : TwObj Φ G) :
    ∃ (B : TwObj Φ G) (φ : A ⟶ B),
      IsIsotropicHull (twPreFrobenioid (G := G) hD hpd) φ :=
  ⟨A, 𝟙 A, (twPreFrobenioid (G := G) hD hpd).Div_id A,
    isPreStep_id _ A, twIsotropic hD hpd A,
    fun Cc _ γ => ⟨γ, (Category.id_comp γ).symm, fun β hβ => by
      have hg : γ = β := by simpa using hβ
      exact hg.symm⟩⟩

/-- ★**(i)(b)** 底の同型は pre-step の span に持ち上がる。

★`𝔽_Φ` の証人に `G` 成分 `1` を添えるだけ ——
★**結論の等式は `𝒟` の中の話なので `G` 成分に影響されない。** -/
theorem twPreStepSpan (A B : TwObj Φ G)
    (α : ((twPreFrobenioid (G := G) hD hpd).toElem.obj A).base ⟶
      ((twPreFrobenioid (G := G) hD hpd).toElem.obj B).base) (hα : IsIso α) :
    ∃ (X : TwObj Φ G) (φ : X ⟶ A) (ψ : X ⟶ B)
      (hφ : IsPreStep (twPreFrobenioid (G := G) hD hpd) φ),
      IsPreStep (twPreFrobenioid (G := G) hD hpd) ψ ∧
      α = @inv _ _ _ _ ((twPreFrobenioid (G := G) hD hpd).Base φ) hφ.2
            ≫ (twPreFrobenioid (G := G) hD hpd).Base ψ := by
  obtain ⟨X₀, φ₀, ψ₀, hφ₀, hψ₀, heq⟩ :=
    elemFrob_preStepSpan Φ hD hpd A.ofElem B.ofElem α hα
  exact ⟨⟨X₀⟩, ⟨φ₀, 1⟩, ⟨ψ₀, 1⟩, hφ₀, hψ₀, heq⟩

/-- ★★**pull-back 性は `𝔽_Φ` 成分だけで決まる**(次数 1 のとき)。

★**単射性**: `G` 成分は次数 1 なので消約される。
★**全射性**: 目標の `G` 成分から `f.unit⁻¹` を掛けて逆算する。 -/
theorem twIsPullBack_of {A B : TwObj Φ G} (f : A ⟶ B) (hlin : f.hom.deg = 1)
    (h : IsPullBack (elemPreFrobenioid Φ hD hpd) f.hom) :
    IsPullBack (twPreFrobenioid (G := G) hD hpd) f := by
  intro X
  obtain ⟨hinj, hsurj⟩ := h X.ofElem
  constructor
  · intro g₁ g₂ hg
    have hv := congrArg Subtype.val hg
    have h1 : g₁ ≫ f = g₂ ≫ f := congrArg Prod.fst hv
    have h2 : g₁.hom.base = g₂.hom.base := congrArg Prod.snd hv
    have hhom : g₁.hom = g₂.hom :=
      hinj (Subtype.ext (Prod.ext (congrArg TwHom.hom h1) h2))
    refine TwHom.ext hhom ?_
    have hu := congrArg TwHom.unit h1
    rw [twComp_unit, twComp_unit, hlin] at hu
    simpa using mul_right_cancel hu
  · intro y
    obtain ⟨f₀, hf₀⟩ := hsurj ⟨(y.1.1.hom, y.1.2), y.2⟩
    have hf₀' := congrArg Subtype.val hf₀
    have e1 : f₀ ≫ f.hom = y.1.1.hom := congrArg Prod.fst hf₀'
    have e2 : f₀.base = y.1.2 := congrArg Prod.snd hf₀'
    refine ⟨⟨f₀, y.1.1.unit * f.unit⁻¹⟩, Subtype.ext (Prod.ext ?_ e2)⟩
    refine TwHom.ext e1 ?_
    show (y.1.1.unit * f.unit⁻¹) ^ ((f.hom.deg : ℕ+) : ℕ) * f.unit = y.1.1.unit
    rw [hlin]
    simp

/-- ★**(iv)(b)** pull-back は LB-invertible かつ linear。

★co-angular は無条件、isometric と linear は `𝔽_Φ` 側から。 -/
theorem twPullBackLB {A B : TwObj Φ G} (α : A ⟶ B)
    (h : IsPullBack (twPreFrobenioid (G := G) hD hpd) α) :
    IsLBInvertible (twPreFrobenioid (G := G) hD hpd) α ∧
      IsLinear (twPreFrobenioid (G := G) hD hpd) α := by
  -- `𝔽_Φ` 側の pull-back 性を取り出す
  have h0 : IsPullBack (elemPreFrobenioid Φ hD hpd) α.hom := by
    intro X₀
    obtain ⟨hinj, hsurj⟩ := h ⟨X₀⟩
    constructor
    · intro g₁ g₂ hg
      have hv := congrArg Subtype.val hg
      have e1 : g₁ ≫ α.hom = g₂ ≫ α.hom := congrArg Prod.fst hv
      have e2 : g₁.base = g₂.base := congrArg Prod.snd hv
      have := hinj (a₁ := (⟨g₁, 1⟩ : (⟨X₀⟩ : TwObj Φ G) ⟶ A))
        (a₂ := (⟨g₂, 1⟩ : (⟨X₀⟩ : TwObj Φ G) ⟶ A))
        (Subtype.ext (Prod.ext (TwHom.ext e1 rfl) e2))
      exact congrArg TwHom.hom this
    · intro y
      obtain ⟨g, hg⟩ := hsurj ⟨((⟨y.1.1, 1⟩ : (⟨X₀⟩ : TwObj Φ G) ⟶ B), y.1.2), y.2⟩
      have hg' := congrArg Subtype.val hg
      have hp1 : g ≫ α = (⟨y.1.1, 1⟩ : (⟨X₀⟩ : TwObj Φ G) ⟶ B) := congrArg Prod.fst hg'
      have hp2 : g.hom.base = y.1.2 := congrArg Prod.snd hg'
      exact ⟨g.hom, Subtype.ext (Prod.ext (congrArg TwHom.hom hp1) hp2)⟩
  obtain ⟨hlb, hlin⟩ := elemFrob_pullBackLB Φ hD hpd α.hom h0
  exact ⟨⟨twCoAngular hD hpd α, hlb.2⟩, hlin⟩

end Bridge

/-! ### ★測定 —— ここまでで確定したこと(2026-08-16)

★**確定した**:
- 捻れ積は**圏である**(結合律・単位律を証明した)
- 射影関手 `twProj : 𝔽_Φ ⋉ G ⥤ 𝔽_Φ` がある
- ★★**次数 `d` の射(`G` 成分が `1`)は、`G` に `d`-捻れがあれば mono でない**
- ★**totally epimorphic である**(`G` 成分は群の消約)
- ★**connected である**(`𝔽_Φ` の証明をそのまま写せる)
- ★★**`PreFrobenioid` 構造を持つ**(`twPreFrobenioid`)

★**未確定(残る仕事)**:
- ★**`Definition 1.3`(`FrobenioidCore` ＋ `Frobenioid`)の全条件を満たすか** ——
  これが済めば `Gap_1_14_iii` は ③(`sourceGap`)に上がる。
  ★`elemFrob_isFrobenioid`(`Proposition 1.5`)の各条件に `G` 成分の処理を
  足す形になる
- `Frobenius-normalized` が実際に成り立つこと(docstring で計算したが未形式化)

★★**ここまでで「実物の圏で mono が破れる」ことは確定した。**
★**残るのは「その圏が Frobenioid である」ことだけである。**

★★**`sorry` は置かない。** 未確定のものは**書かない**。
-/

end ABC3.Check.FrdI
