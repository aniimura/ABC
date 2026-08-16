import ABC3.Found.GenEll.LogDiffValue
import Mathlib.AlgebraicGeometry.IdealSheaf.Functorial

/-!
# [GenEll] Definition 1.5, (iv) —— log-cond(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> determines a well-deﬁned log-diﬀerent function log-diﬀX on X(Q).

## ★★自分の測定を訂正する

本日(2026-08-17)私は `blocked-leaves.json` に

> 「`grep -rln Cartier Mathlib/` が **0 件**。……『重い』のではなく『**底が抜けている**』」

と書いた。★**この判断は誤りだった。**

`Cartier` という**名前**は確かに mathlib に 1 件も無い。
しかし原文が要求している**中身**——「因子 `D ⊆ X` を `x_F : Spec(𝓞_F) → X` で引き戻す」——は

- **`AlgebraicGeometry.Scheme.IdealSheafData.comap`**(イデアル層の引き戻し)
- **`Scheme.IdealSheafData.equivOfIsAffine`**(アフィンならイデアル層 ≃ イデアル)
- **`Scheme.ΓSpecIso`**(`Γ(Spec R, ⊤) ≅ R`)

の 3 つで**そのまま届く**。実際、本ファイルはそれだけで `f_x^D` を構成する。

★★**名前で測って「無い」と結論したのが誤りの型である。**
★これで「正面から要ると思ったものが要らなかった」の 6 例目になる。

## ★何を取るか

- `IsEffectiveCartier` —— 有効 Cartier 因子(イデアル層が局所的に**非零因子 1 つ**で生成)
- `pullbackIdeal` —— `D_x`(`D` を `x_F` で引き戻したイデアル)
- `conductorADiv` —— `f_x^D ≝ (D_x)_red`(原文の被約化は**根基**である)
- `logCond` —— `log-cond_D(x) ≝ deg_F(f_x^D)`

## ★負の対照

`D = ⊤`(空因子)なら `log-cond = 0` でなければならない(`logCond_top`)。
★これが出ないなら引き戻しの向きか被約化を取り違えている。

## ★★まだ届いていないもの

`x_F` を「点 `x ∈ X(F)` から作る」ところ(`Definition 1.5, (i)` の
minimal field of definition)は本ファイルの外である。
★本ファイルは **`x_F` を与えられたものとして受ける**——
そこは `Scheme.Hom.image` があるので道はあるが、まだ歩いていない。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField IsDedekindDomain

/-! ## ★有効 Cartier 因子 -/

/-- **有効 Cartier 因子**。

イデアル層 `I` が、各点のまわりのあるアフィン開集合の上で
**非零因子 1 つ**で生成されること(Stacks 01WQ の定義)。

★`Cartier` は mathlib に無いので、ここで定義する。
★「局所的に主イデアル」だけでは足りない——**非零因子**でなければ
`D` が余次元 1 にならない。 -/
def IsEffectiveCartier {X : Scheme} (I : X.IdealSheafData) : Prop :=
  ∀ x : X, ∃ U : X.affineOpens, x ∈ U.1 ∧ ∃ f : Γ(X, U),
    I.ideal U = Ideal.span {f} ∧ f ∈ nonZeroDivisors Γ(X, U)

/-- ★**空因子は有効 Cartier である**(非空虚性の witness)。

`⊤ : IdealSheafData`(= 空の閉部分スキーム)は `1` で生成され、`1` は非零因子。
★これが無いと `IsEffectiveCartier` が空虚な述語でないことを誰も保証しない。 -/
theorem isEffectiveCartier_top (X : Scheme) :
    IsEffectiveCartier (⊤ : X.IdealSheafData) := by
  intro x
  obtain ⟨_, ⟨V, hV, rfl⟩, hxV, -⟩ :=
    X.isBasis_affineOpens.exists_subset_of_mem_open (Set.mem_univ x) isOpen_univ
  refine ⟨⟨V, hV⟩, hxV, 1, ?_, ?_⟩
  · rw [Ideal.span_singleton_one]
    rfl
  · exact Submonoid.one_mem _

/-! ## ★引き戻し -/

section Pullback

variable (F : Type*) [Field F] [NumberField F]

/-- `Spec(𝓞_F)` をスキームとして。原文の `Spec(OF )`。 -/
noncomputable abbrev specRingOfIntegers : Scheme := Spec (CommRingCat.of (𝓞 F))

variable {X : Scheme}

/-- ★★**`D_x`** —— 因子 `D ⊆ X` を `x_F : Spec(𝓞_F) → X` で引き戻したイデアル。

原文 (GenEll p.8):
> Then the morphism Spec(OF ) → X determined by x [where we recall that X is proper!] allows one to pull-back the divisor D to Spec(OF ) so as to obtain an eﬀective divisor Dx on Spec(OF ).

★3 段の合成である:
`IdealSheafData.comap`(引き戻し)→ `equivOfIsAffine`(アフィンなのでイデアル)
→ `ΓSpecIso`(`Γ(Spec 𝓞_F, ⊤) ≅ 𝓞_F`)。 -/
noncomputable def pullbackIdeal (D : X.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) : Ideal (𝓞 F) :=
  (Scheme.IdealSheafData.equivOfIsAffine (D.comap xF)).comap
    (Scheme.ΓSpecIso (CommRingCat.of (𝓞 F))).inv.hom

/-- ★**空因子の引き戻しは単位イデアル**。

`Scheme.IdealSheafData.comap_top` から出る。 -/
@[simp] theorem pullbackIdeal_top (xF : specRingOfIntegers F ⟶ X) :
    pullbackIdeal F (⊤ : X.IdealSheafData) xF = ⊤ := by
  simp [pullbackIdeal, Scheme.IdealSheafData.comap_top]

/-! ## ★★導手と log-cond -/

/-- ★★**`f_x^D ≝ (D_x)_red`** —— 導手。

原文 (GenEll p.8):
> Note that Dx, (Dx)red may also be regarded as arithmetic divisors ∈ ADiv(F ) that are supported in V(F )non.

★被約化 `(−)_red` は、イデアルの側では**根基**である
(`Spec(𝓞_F)` の閉部分スキームの被約構造 ↔ 根基イデアル)。 -/
noncomputable def conductorADiv (D : X.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) : ADiv F :=
  idealADiv F (pullbackIdeal F D xF).radical

/-- ★導手は有効算術因子である(原文 "eﬀective")。 -/
theorem conductorADiv_isEffective (D : X.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) :
    (conductorADiv F D xF).IsEffective :=
  idealADiv_isEffective F _

/-- ★★**`log-cond_D(x) ≝ deg_F(f_x^D)`**(正規化次数)。 -/
noncomputable def logCond (D : X.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) : ℝ :=
  degNormalized (conductorADiv F D xF)

/-- ★★**負の対照** —— `D = ⊤`(空因子)なら `log-cond = 0`。

★引き戻しの向きか被約化を取り違えていれば、ここが 0 にならない。 -/
@[simp] theorem logCond_top (xF : specRingOfIntegers F ⟶ X) :
    logCond F (⊤ : X.IdealSheafData) xF = 0 := by
  have hne : (⊤ : Ideal (𝓞 F)) ≠ 0 := by
    simp only [ne_eq, Ideal.zero_eq_bot]
    exact top_ne_bot
  have : deg (conductorADiv F (⊤ : X.IdealSheafData) xF) = 0 := by
    rw [conductorADiv, pullbackIdeal_top, Ideal.radical_top, deg_idealADiv F _ hne,
      Ideal.absNorm_top]
    simp
  rw [logCond, degNormalized, this, zero_div]

/-- ★★**[GenEll] Proposition 1.6 の非アルキメデス側(引き戻しの言葉で)**。

> `log-cond_D(x) ≤ deg_F(D_x)`

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

★被約化で次数は減る(`deg_idealADiv_radical_le`)。
★★**原文の右辺は「高さ」だが、その手前の段——「導手 ≤ 引き戻した因子の次数」——が
これである。** 高さとの接続には `X` 上の算術直線束(`Definition 1.1`、
複素解析空間が要るので未着手)が要る。**混同しない。** -/
theorem logCond_le_degNormalized_pullback (D : X.IdealSheafData)
    (xF : specRingOfIntegers F ⟶ X) (h : pullbackIdeal F D xF ≠ 0) :
    logCond F D xF ≤ degNormalized (idealADiv F (pullbackIdeal F D xF)) := by
  have hpos : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  rw [logCond, conductorADiv, degNormalized, degNormalized]
  exact div_le_div_of_nonneg_right (deg_idealADiv_radical_le F _ h) hpos.le

end Pullback

/-! ## ★出典の紐付け(`.src`) -/

def logCond.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5, (iv)",
    sectionId := "genell-def-1-5" }

def IsEffectiveCartier.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5, (iv)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
