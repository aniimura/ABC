import ABC3.Found.FrdI.Prop110

/-!
# [FrdI] Proposition 1.11 —— Pull-back and Linear Morphisms

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.36–37(**400 dpi 目視確認 2026-08-16**、親が p.36 / p.37 を描画して照合。
さらに (iii) の 1 行は **600 dpi で拡大確認**した ——下記)。

原文 (FrdI p.36):
> (Pull-back and Linear Morphisms) Let Φ be a divi-

## ★この命題の規模(着手前の測定)

**7 条 (i)–(vii)、主張は 15**:

| 条 | 主張 | 内容 |
|---|---|---|
| (i) | 1 | `𝒞^pl-bk → 𝒟` が **full**(Aut-ample ＋ base-trivial 型) |
| (ii) | 1 | 同じ関手が **faithful**(unit-trivial 型) |
| (iii) | 1 | pull-back に沿った自己射の一意な持ち上げ |
| (iv) | 2 | `𝒪^▷(A) ↪ 𝒪^▷(B)` の存在 ＋ 一意性 |
| (v) | 4 | `Definition 1.3, (iii), (d)` の圏同値の**関手性** |
| (vi) | 4 | pull-back が FSM / fiberwise-surj / mono / irreducible ⟺ その底が |
| (vii) | 2 | co-angular pre-step の持ち上げ ＋「In particular」FSM |

## ★★`Proposition 1.10` の投資が回収されている

- **両方の圏同値**(`coaPreUnderEquiv` / `coaPreOverEquiv`)—— 1.10 で**初めて消費した**
- **irreducible**(§0)—— 1.10 の (iv) のために**このセッションで足した語彙**
- `𝒪^▷` の扱い —— 1.10 の (iii) moreover で慣れた

★**不足は `𝒞^lin`(linear 射のなす広い部分圏)だけ**((v) が要る)。

## ★★★原典の誤植を 1 つ見つけた(600 dpi 目視確認)

原文 (FrdI p.36) の (iii) の末尾は
```
Base(β) = β_𝒟,  α ◦ φ = β ◦ φ.
```
★**`β ◦ φ` は型が通らない。** `β ∈ End_𝒞(B)`、`φ : B → A` なので、
`β ◦ φ`(＝「先に `φ`、次に `β`」)は `φ` の終域 `A` と `β` の始域 `B` が合わない。

★**意図は `α ◦ φ = φ ◦ β`** である。根拠:
- **型が通るのはこれだけ**
- ★**同じ命題の (iv) が `α ◦ φ = φ ◦ β` と書いている**(同じ形の主張)

★★**これは省略でも圧縮でもなく、誤植である。**
我々がこれまで数えてきた 6 種類のギャップとは別の、★**7 種類目: 誤植**。
★**形式化が最初に見つけた「原典の誤り」**である(これまでは全て省略・圧縮だった)。
★**数学的な内容は変わらない** —— 意図は文脈から一意に定まる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(iii) —— pull-back に沿った自己射の一意な持ち上げ

原文 (FrdI p.36):
> (iii) Let φ : B →A be a pull-back morphism that projects to a morphism

原文 (FrdI p.36):
> such that Base(α) ◦φD = φD ◦βD, there exists a unique β ∈EndC(B) such that

★**`IsPullBack` の定義がそのまま (iii) である**:
```
IsPullBack φ := ∀ X, Bijective (f ↦ (f ≫ φ, Base f))
```
★**全射性が存在を、単射性が一意性を与える。**
★**原文が `Definition 1.2, (ii)` として定義したものを、(iii) は
「自己射の持ち上げ」の形で述べ直しているだけ**である。
-/

include P in
/-- **(iii)** —— pull-back `φ : B ⟶ A` に沿って、底の自己射が一意に持ち上がる。

★**原典の誤植を訂正した形で述べる**(上のファイル docstring を見よ):
原文は `α ◦ φ = β ◦ φ` と書くが型が通らず、意図は `α ◦ φ = φ ◦ β`
(＝ Lean の `φ ≫ α = β ≫ φ`)である。 -/
theorem prop_1_11_iii {A B : C} (φ : B ⟶ A) (hφ : IsPullBack P φ)
    (α : End A) (βD : End ((P.toElem.obj B).base))
    (hsq : P.Base φ ≫ P.Base (α : A ⟶ A) = (βD : _ ⟶ _) ≫ P.Base φ) :
    ∃! β : End B, P.Base (β : B ⟶ B) = βD ∧
      (φ ≫ (α : A ⟶ A)) = (β : B ⟶ B) ≫ φ := by
  obtain ⟨hinj, hsurj⟩ := hφ B
  obtain ⟨β, hβ⟩ := hsurj ⟨(φ ≫ (α : A ⟶ A), (βD : _ ⟶ _)), by
    rw [P.Base_comp, hsq]⟩
  have hβ' := congrArg Subtype.val hβ
  refine ⟨β, ⟨?_, ?_⟩, ?_⟩
  · exact congrArg Prod.snd hβ'
  · exact (congrArg Prod.fst hβ').symm
  · rintro γ ⟨hγ1, hγ2⟩
    have hf : (γ : B ⟶ B) ≫ φ = (β : B ⟶ B) ≫ φ := by
      rw [← hγ2]
      exact (congrArg Prod.fst hβ').symm
    have hb : P.Base (γ : B ⟶ B) = P.Base (β : B ⟶ B) := by
      rw [hγ1]
      exact (congrArg Prod.snd hβ').symm
    exact hinj (Subtype.ext (Prod.ext hf hb))

/-! ## ★(ii) —— `𝒞^pl-bk → 𝒟` の faithful 性

原文 (FrdI p.36):
> (ii) Suppose further that C is of unit-trivial type. Then the natural projection

★**原文の証明**(p.36–37、目視)は 4 段:
1. `φ, ψ : A ⟶ B` が同じ底へ落ちる pull-back 射なら、
   **pull-back の定義から** base-identity 自己射 `α, β ∈ End(A)` が
   `ψ = φ ∘ α`、`φ = ψ ∘ β` を満たす
2. すると `ψ = ψ ∘ β ∘ α`、`φ = φ ∘ α ∘ β` なので、
   **再び pull-back の定義(の一意性)から** `α ∘ β = β ∘ α = 𝟙`
3. よって `α, β ∈ Aut_𝒞(A)`、したがって `α, β ∈ 𝒪^×(A)`
4. **unit-trivial 型**なので `𝒪^×(A) = {1}`、よって `φ = ψ`

★★**`Definition 1.2, (ii)` を 2 回使う**(存在と一意性)のが要点である。
★**原文は「it thus follows formally」と書くが、2 回使うことは書いていない。**
-/

include P in
/-- **(ii)** —— unit-trivial 型なら、同じ底へ落ちる pull-back 射は一致する
(＝ `𝒞^pl-bk → 𝒟` が faithful)。 -/
theorem prop_1_11_ii {A B : C} (hut : IsUnitTrivial P A)
    (φ ψ : A ⟶ B) (hφ : IsPullBack P φ) (hψ : IsPullBack P ψ)
    (hbase : P.Base φ = P.Base ψ) : φ = ψ := by
  -- 段1: `α` と `β` を取る
  obtain ⟨hφinj, hφsurj⟩ := hφ A
  obtain ⟨hψinj, hψsurj⟩ := hψ A
  obtain ⟨α, hα⟩ := hφsurj ⟨(ψ, 𝟙 _), by rw [hbase, Category.id_comp]⟩
  obtain ⟨β, hβ⟩ := hψsurj ⟨(φ, 𝟙 _), by rw [← hbase, Category.id_comp]⟩
  have hα' := congrArg Subtype.val hα
  have hβ' := congrArg Subtype.val hβ
  have hαφ : α ≫ φ = ψ := congrArg Prod.fst hα'
  have hαb : P.Base α = 𝟙 _ := congrArg Prod.snd hα'
  have hβψ : β ≫ ψ = φ := congrArg Prod.fst hβ'
  have hβb : P.Base β = 𝟙 _ := congrArg Prod.snd hβ'
  -- 段2: `α ≫ β = 𝟙`(★pull-back の**一意性**を使う)
  have hcomp : β ≫ α = 𝟙 A := by
    refine hφinj (Subtype.ext (Prod.ext ?_ ?_))
    · show (β ≫ α) ≫ φ = 𝟙 A ≫ φ
      rw [Category.assoc, hαφ, hβψ, Category.id_comp]
    · show P.Base (β ≫ α) = P.Base (𝟙 A)
      rw [P.Base_comp, hαb, hβb, Category.id_comp, P.Base_id]
  have hcomp' : α ≫ β = 𝟙 A := by
    refine hψinj (Subtype.ext (Prod.ext ?_ ?_))
    · show (α ≫ β) ≫ ψ = 𝟙 A ≫ ψ
      rw [Category.assoc, hβψ, hαφ, Category.id_comp]
    · show P.Base (α ≫ β) = P.Base (𝟙 A)
      rw [P.Base_comp, hβb, hαb, Category.id_comp, P.Base_id]
  -- 段3: `α ∈ 𝒪^×(A)`
  haveI : IsIso α := ⟨β, hcomp', hcomp⟩
  have hmem : (α : End A) ∈ OTimes P A := by
    refine ⟨⟨?_, ?_⟩, (CategoryTheory.isUnit_iff_isIso (α : End A)).mpr inferInstance⟩
    · show P.Base α = P.Base (𝟙 A)
      rw [hαb, P.Base_id]
    · exact degFr_of_isIso P α
  -- 段4: unit-trivial から `α = 𝟙`
  rw [hut] at hmem
  have hone : α = 𝟙 A := hmem
  rw [← hαφ, hone, Category.id_comp]

/-! ## ★(iv) —— `𝒪^▷(A) ↪ 𝒪^▷(B)`

原文 (FrdI p.36):
> (iv) Every co-angular linear morphism φ : B →A determines an injection

★★**逐語の注意(事故 #5 の再確認)**: PDF は `𝒪^▷(A) ↪ 𝒪^▷(B)`(**単射の矢印**)だが
`pdftotext` は `→` に変える。★**400 dpi 目視で確認した。**
散文に「an injection」とあるので意味は復元できるが、**矢印だけを見ると単射性が消える。**

★**原文の証明**(p.37、目視)は 3 つに分かれる:
> the existence of the map … follows immediately from assertion (iii) and Definition 1.3, (iii), (c);
> **the asserted injectivity of this map follows from the total epimorphicity of C**;
> the fact that this map is uniquely determined … follows from the fact that
> **pre-steps are monomorphisms** … and the definition of a "pull-back morphism"

★★**存在・単射性・一意性が 3 つの別々の理由から来る。**
★**そして単射性の理由は `𝒞` の totally epimorphicity である** ——
★`Proposition 1.10, (iii)` で我々が見つけた **epi による消去**が、ここでも要点になっている。
-/

include P in
/-- ★**(iv) の単射性** —— 対応 `α ↦ β` は単射。

★**理由は `𝒞` の totally epimorphicity**(原文どおり)。
`φ ≫ α₁ = β ≫ φ = φ ≫ α₂` から `φ` が epi で `α₁ = α₂`。 -/
theorem prop_1_11_iv_injective {A B : C} (φ : B ⟶ A)
    {α₁ α₂ : End A} {β : End B}
    (h₁ : φ ≫ (α₁ : A ⟶ A) = (β : B ⟶ B) ≫ φ)
    (h₂ : φ ≫ (α₂ : A ⟶ A) = (β : B ⟶ B) ≫ φ) : α₁ = α₂ := by
  haveI : Epi φ := P.totEpiC _ _ _
  exact (cancel_epi φ).mp (h₁.trans h₂.symm)

include P in
/-- ★**(iv) の一意性** —— `α` に対応する `β` は高々 1 つ。

★**理由は pre-step が mono であること**(原文どおり、`Definition 1.3, (v), (a)`)。
★**単射性とは違う理由から来る** —— 単射性は epi、一意性は mono。
★★**同じ命題の中で、消去の向きが 2 回とも別々に要る。** -/
theorem prop_1_11_iv_unique (F : FrobenioidCore P) {A B : C} (φ : B ⟶ A)
    (hφs : IsPreStep P φ) {α : End A} {β₁ β₂ : End B}
    (h₁ : φ ≫ (α : A ⟶ A) = (β₁ : B ⟶ B) ≫ φ)
    (h₂ : φ ≫ (α : A ⟶ A) = (β₂ : B ⟶ B) ≫ φ) : β₁ = β₂ := by
  haveI : Mono φ := F.preStepMono φ hφs
  exact (cancel_mono φ).mp (h₁.symm.trans h₂)

include P in
/-- ★★**(iv) の存在** —— co-angular linear `φ : B ⟶ A` に沿って
`𝒪^▷(A)` の元が `𝒪^▷(B)` へ移る。

★**構成は原文どおり 3 手**:
1. `Proposition 1.7, (iii)` で `φ = β ≫ α`(`β` pre-step、`α` pull-back)と分解
   ——★**`β` が co-angular であることは `Proposition 1.4, (iv)` から**
2. ★**(iii)** を pull-back `α` に当てて `γ_X ∈ End X` を得る
   (★**`γ_X` が linear であることは次数の消去から出る** ——原文は書いていない)
3. ★**`otriBwd`**(`Definition 1.3, (iii), (c)`)を co-angular pre-step `β` に当てて
   `δ ∈ 𝒪^▷(B)` を得る -/
theorem prop_1_11_iv_exists (F : FrobenioidCore P) {A B : C} (φ : B ⟶ A)
    (hφc : IsCoAngular P φ) (hφl : IsLinear P φ) (γA : End A) (hγA : γA ∈ OTri P A) :
    ∃ δ : End B, δ ∈ OTri P B ∧ φ ≫ (γA : A ⟶ A) = (δ : B ⟶ B) ≫ φ := by
  -- 手1: 分解
  obtain ⟨X, β, α, hfac, hβs, hαpb⟩ := (prop_1_7_iii_linear_factor P F φ).mp hφl
  have hβc : IsCoAngular P β :=
    prop_1_4_iv_mpr P F (show φ = 𝟙 B ≫ β ≫ α by rw [hfac, Category.id_comp])
      (isFrobeniusType_of_isIso P (𝟙 B)) hβs hαpb hφc
  -- 手2: (iii) を `α` に当てる
  have hsq : P.Base α ≫ P.Base (γA : A ⟶ A) = 𝟙 _ ≫ P.Base α := by
    rw [show P.Base (γA : A ⟶ A) = P.Base (𝟙 A) from hγA.1, P.Base_id,
      Category.comp_id, Category.id_comp]
  obtain ⟨γX, ⟨hγXb, hγXsq⟩, _⟩ := prop_1_11_iii P α hαpb γA (𝟙 _) hsq
  -- `γ_X` は linear(★次数の消去。原文は書いていない)
  have hγXl : IsLinear P (γX : X ⟶ X) := by
    have h := congrArg P.degFr hγXsq
    rw [P.degFr_comp, P.degFr_comp, show P.degFr (γA : A ⟶ A) = 1 from hγA.2,
      one_mul] at h
    show P.degFr (γX : X ⟶ X) = 1
    have h2 : P.degFr α * 1 = P.degFr α * P.degFr (γX : X ⟶ X) := by
      rw [mul_one]; exact h
    exact (mul_left_cancel h2).symm
  have hγXm : (γX : End X) ∈ OTri P X := by
    refine ⟨?_, hγXl⟩
    show P.Base (γX : X ⟶ X) = P.Base (𝟙 X)
    rw [hγXb, P.Base_id]
  -- 手3: `otriBwd`
  obtain ⟨δ, ⟨hδm, hδsq⟩, _⟩ := F.otriBwd β hβc hβs (γX : End X) hγXm
  refine ⟨δ, hδm, ?_⟩
  rw [hfac, Category.assoc, hγXsq, ← Category.assoc, hδsq, Category.assoc]

/-! ## ★(i) —— `𝒞^pl-bk → 𝒟` の full 性

原文 (FrdI p.36):
> (i) Suppose further that C is of Aut-ample and base-trivial type. Then the

★**原文の証明**(p.37、目視)は 4 手:
1. `Definition 1.3, (i), (c)` の圏同値の**本質的全射性**で、
   `φ_𝒟` に対応する pull-back 射 `ψ : C ⟶ B` を得る(`C_𝒟 ≅ A_𝒟`)
2. ★**base-trivial 型**なので `C ≅ A`、よって `A = C` としてよい
3. すると `φ_𝒟 = ψ_𝒟 ∘ δ`(`δ ∈ Aut_𝒟(A_𝒟)`)と書ける
4. ★**Aut-ample 型**なので `δ` が `γ ∈ Aut_𝒞(A)` に持ち上がり、`ψ ∘ γ` が求めるもの

★★**2 つの型が別々の手で効く** —— base-trivial は**対象**を合わせ、
Aut-ample は**射**を合わせる。★**構造化のときの見立てが当たった。**

★**ここでは手 2–4 を実装する**(手 1 の本質的全射性は `plBkOverFunctor` の
プランビングで、`coaPre_realize` と同型の作業)。★**手 1 の出力を仮定として型に出す。**
-/

include P in
/-- ★**(i) の手 2–4** —— pull-back `ψ : Cc ⟶ B` と底の同型が与えられたとき、
base-trivial ＋ Aut-ample から `A ⟶ B` の pull-back 射で底が `φ_𝒟` のものを作る。

★**仮定 `hψ` が手 1(圏同値の本質的全射性)の出力**である。 -/
theorem prop_1_11_i_lift (F : FrobenioidCore P) {A B Cc : C}
    (hbt : IsBaseTrivial P A) (haa : IsAutAmple P A)
    (φD : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base)
    (ψ : Cc ⟶ B) (hψpb : IsPullBack P ψ)
    (θ : (P.toElem.obj Cc).base ≅ (P.toElem.obj A).base)
    (hθ : P.Base ψ = θ.hom ≫ φD) :
    ∃ χ : A ⟶ B, IsPullBack P χ ∧ P.Base χ = φD := by
  -- 手2: base-trivial で `Cc ≅ A`
  obtain ⟨ι⟩ := hbt Cc ⟨θ.symm⟩
  -- 手3: `δ := Base ι.inv ≫ θ.hom` は `Aut_𝒟(A_𝒟)`
  haveI : IsIso (P.Base ι.inv) := by
    refine ⟨P.Base ι.hom, ?_, ?_⟩
    · rw [← P.Base_comp, ι.inv_hom_id, P.Base_id]
    · rw [← P.Base_comp, ι.hom_inv_id, P.Base_id]
  set δ : End ((P.toElem.obj A).base) := P.Base ι.inv ≫ θ.hom with hδ
  haveI : IsIso (δ : _ ⟶ _) := by rw [hδ]; infer_instance
  -- 手4: Aut-ample で `δ⁻¹` を持ち上げる
  obtain ⟨γ, hγiso, hγb⟩ := haa (inv (δ : _ ⟶ _)) inferInstance
  haveI : IsIso (γ : A ⟶ A) := hγiso
  refine ⟨(γ : A ⟶ A) ≫ ι.inv ≫ ψ, ?_, ?_⟩
  · exact IsPullBack.comp P (isPullBack_of_isIso P _)
      (IsPullBack.comp P (isPullBack_of_isIso P _) hψpb)
  · rw [P.Base_comp, P.Base_comp, hγb, hθ]
    calc inv (δ : _ ⟶ _) ≫ P.Base ι.inv ≫ θ.hom ≫ φD
        = (inv (δ : _ ⟶ _) ≫ (δ : _ ⟶ _)) ≫ φD := by
          simp only [hδ, Category.assoc]
      _ = φD := by rw [IsIso.inv_hom_id, Category.id_comp]

/-! ## ★(vi) —— pull-back の 4 性質は底と対応する

原文 (FrdI p.37):
> (vi) A pull-back morphism φ ∈Arr(C) is an FSM-morphism (respectively,

原文 (FrdI p.37):
> and only if Base(φ) ∈Arr(D) is.

★**原文の証明**(p.38、目視):
> assertion (vi) follows **formally** from the isomorphism of functors appearing in the
> definition of a "pull-back morphism" [cf. Definition 1.2, (ii)], together with the
> **equivalence of categories induced by the projection functor in Definition 1.3, (i), (c)**

★★**着手前の予測「(vi) は要注意」が当たった。** 4 つの性質のうち、
★**`Definition 1.2, (ii)` だけで出るのは `mono` の `⟸` 向きだけ**である。
他の 3 つと `mono` の `⟹` 向きは、★**`𝒟` の対象・射を `𝒞` へ持ち上げる**必要があり、
それには `Definition 1.3, (i), (c)` の圏同値(＝ (i) の full 性)が要る。

★★**「formally」はここでも「新しい材料が要らない」の意味である** ——
`Definition 1.2, (ii)` と `Definition 1.3, (i), (c)` の 2 つで足りる。
★**しかし「どちらがどの向きに要るか」は書かれていない。**

★**ここでは `mono` の `⟸` 向きを実装する**(`Definition 1.2, (ii)` だけで出る)。
-/

include P in
/-- ★**(vi) の mono の `⟸` 向き** —— 底が mono なら pull-back 射も mono。

★★**`Definition 1.2, (ii)` の単射性だけで出る。** 圏同値は要らない。
★**`Proposition 1.11, (i)` の full 性が要るのは残りの向きである。** -/
theorem prop_1_11_vi_mono_of_baseMono {A B : C} (φ : A ⟶ B) (hφ : IsPullBack P φ)
    (hm : Mono (P.Base φ)) : Mono φ := by
  refine ⟨fun {X} f g h => ?_⟩
  obtain ⟨hinj, _⟩ := hφ X
  refine hinj (Subtype.ext (Prod.ext h ?_))
  show P.Base f = P.Base g
  refine (cancel_mono (P.Base φ)).mp ?_
  rw [← P.Base_comp, ← P.Base_comp, h]

/-! ## ★(vii) —— co-angular pre-step は FSM 射

原文 (FrdI p.37):
> (vii) Let φ : A →B be a co-angular pre-step;

★**続きは引用できない**: `ϵ`(varepsilon)が layout モードで消える。
400 dpi 目視では「`ϵ : C → B` を任意の射とする」と続く。
★**構造化 HTML では `data-txt` で明示したが、Lean の docstring には
その逃げ道が無い**(第27段で記録した制約)。

★★**主張本体は `IsFiberwiseSurjective φ` そのものである。**
`IsFiberwiseSurjective (β : B ⟶ A) := ∀ Z (γ : Z ⟶ A), ∃ D δB δZ, δB ≫ β = δZ ≫ γ`
と (vii) の「`ϵ ∘ γ = φ ∘ α` を満たす `γ`・`α` が存在する」は**同じ形**である
(★ただし (vii) は `γ` が **co-angular pre-step** であることまで言うので**強い**)。

★**原文の証明**(p.38、目視)は `ϵ` について **4 場合に分ける**
(pull-back / isometric pre-step / co-angular pre-step / Frobenius 型)。
分け方は `Definition 1.3, (iv), (a)` と `(v), (b)` の分解による。

★**ここでは「In particular」を実装する** —— 主張本体を仮定に取れば
`IsFSMMorphism = IsFiberwiseSurjective ∧ Mono` の残り半分は
★**`preStepMono`(`Definition 1.3, (v), (a)`)から直ちに出る。**
-/

include P in
/-- ★**(vii) の「In particular」** —— co-angular pre-step は FSM 射。

★**`IsFSMMorphism = IsFiberwiseSurjective ∧ Mono`** なので、
主張本体(fiberwise-surjectivity)が与えられれば、
残りは `preStepMono` から直ちに出る。

★**原文の「In particular」は、定義の分解そのものである** ——
`Proposition 1.10, (iv)` の「In particular」(無限性)とは種類が違う。 -/
theorem prop_1_11_vii_fsm (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hφs : IsPreStep P φ) (hfs : IsFiberwiseSurjective φ) : IsFSMMorphism φ :=
  ⟨hfs, F.preStepMono φ hφs⟩

end ABC3.Found.FrdI
