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
