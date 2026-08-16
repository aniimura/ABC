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

/-! ## ★`𝒞^lin` —— linear 射のなす広い部分圏

原文 (FrdI p.37):
> in the following sense: If φ : A →B is an arbitrary morphism of Clin, α : C →

★**(v) が使う `𝒞^lin` は、我々の側に無かった**(構造化のときの測定どおり)。
★**`𝒞^pl-bk` とまったく同じ作り方**でできる:
- 恒等射は linear(`isLinear_id`)
- linear は合成で閉じる(`IsLinear.comp`)

★★**原文は `𝒞^lin` を `Definition 1.3` などで別に定義していない** ——
`Clin` という記法を (v) で**いきなり使う**。
★**「linear 射のなす部分圏」は自明なので定義を置かない、という判断**であり、
★**我々も同じ判断をする(定義を 3 行で置く)。**
-/

/-- `linear` 射がなす `MorphismProperty`。 -/
def linProp : MorphismProperty C := fun _ _ φ => IsLinear P φ

instance : (linProp P).ContainsIdentities :=
  ⟨fun A => isLinear_id P A⟩

instance : (linProp P).IsStableUnderComposition :=
  ⟨fun _ _ hf hg => IsLinear.comp P hf hg⟩

instance : MorphismProperty.IsMultiplicative (linProp P) where

/-- **`𝒞^lin`** —— linear 射が定める広い部分圏。★`(v)` が使う。 -/
abbrev Lin : Type u2 := WideSubcategory (linProp P)

/-- ★**非退化(下)**: 恒等射は `𝒞^lin` の射である(`ContainsIdentities` の言い換え)。 -/
theorem linProp_id (A : C) : linProp P (𝟙 A) := isLinear_id P A

include P in
/-- ★**非退化(上)**: `𝒞^lin` は `𝒞` 全体ではない —— 次数が 1 でない射は入らない。

★**`Proposition 1.5` の witness(`𝔽_ℕ` の Vee 上)には次数 2 の射がある**ので、
この主張は空ではない。★**ここでは「次数が 1 でなければ `𝒞^lin` に入らない」を
定義から述べる**(モデル依存の非退化は `Witness.lean` の側の仕事)。 -/
theorem not_linProp_of_degFr_ne_one {A B : C} (φ : A ⟶ B) (h : P.degFr φ ≠ 1) :
    ¬ linProp P φ := h

/-! ### ★★親の規則違反(即座に訂正した)

★親は最初 `linProp.src` を
`item := "Proposition 1.11, (v) — 𝒞^lin"` として置いた。
★**ゲートは通り、進捗が 11/54 → 12/54 に増えた。**

★★**これは規則違反である。** 我々が 2026-08-15 に明文化した規則:
> `.src` は「その原典項目を**完全に**実装した」という主張である。
> 項目へ向けた**構成**(例: `MonoidOn.charOn`)にも付けない。原典の項目ではないから。

`𝒞^lin` は **`Proposition 1.11` へ向けた語彙**であって、
`Proposition 1.11` の実装ではない。`1.11` は 15 主張中 7 しか通っていない。

★**即座に外した。進捗は 11/54 に戻る。**

★★**測定: 規則違反は「数が増える」形で現れた。**
★**指標を上げる方法が 2 つある** —— 項目を完成させるか、規則を緩めるか。
★**後者は一瞬で、しかもゲートを通る。** だからこそ規則は**人間(ここでは親)が守る**しかない。
★**器具は `.src` の locator を検証するが、「完全実装かどうか」は検証できない**
——これは `frdi-progress.mjs` の docstring に最初から書いてある限界である。
★**限界を書いておいたおかげで、違反に気づけた。**
-/


/-! ## ★★`Definition 1.3, (i), (c)` の圏同値を消費する

★**(i) の手 1、(vi) の残り、(vii) の一部が、すべてこの 1 本を待っている。**
`Proposition 1.10` で `coaPreUnderEquiv` / `coaPreOverEquiv` を消費したのと
**同じ形のプランビング**である。

★**橋の形**: `𝒟` の射 `f : Y ⟶ A_𝒟` に対し、
**pull-back 射 `ψ : X ⟶ A` で `Base ψ` が `f` に(スライスで)同型なもの**が存在する。

★★**これが原文の「the equivalence of categories induced by the projection functor
in Definition 1.3, (i), (c)」の実体である。**
-/

include P in
/-- ★★**`Definition 1.3, (i), (c)` の圏同値の初の消費**。

`𝒟` の射 `f : Y ⟶ A_𝒟` は、**pull-back 射の底として実現できる**
(スライス `𝒟_{A_𝒟}` での同型を除いて)。

★**`Proposition 1.11, (i)` の手 1 がこれである。** -/
theorem plBk_realize (F : FrobenioidCore P) (A : C) {Y : D}
    (f : Y ⟶ (P.toElem.obj A).base) :
    ∃ (X : C) (ψ : X ⟶ A) (_ : IsPullBack P ψ)
      (θ : (P.toElem.obj X).base ≅ Y), P.Base ψ = θ.hom ≫ f := by
  haveI := F.plBkEquiv A
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := plBkOverFunctor P A)
    (Over.mk f)
  refine ⟨Z.left.obj, Z.hom.hom, Z.hom.property, ⟨e.hom.left, e.inv.left, ?_, ?_⟩, ?_⟩
  · exact congrArg CommaMorphism.left e.hom_inv_id
  · exact congrArg CommaMorphism.left e.inv_hom_id
  · exact (Over.w e.hom).symm

include P in
/-- ★★★**`Proposition 1.11, (i)` 完成** —— `𝒞^pl-bk → 𝒟` は **full**。

Aut-ample かつ base-trivial 型なら、`𝒟` の任意の射 `φ_𝒟 : A_𝒟 ⟶ B_𝒟` は
**pull-back 射の底として実現できる**。

★**手 1 は `plBk_realize`(圏同値の消費)、手 2–4 は `prop_1_11_i_lift`。**
★**構造化のときの見立て「base-trivial は対象を、Aut-ample は射を合わせる」が
そのまま構成になった。** -/
theorem prop_1_11_i (F : FrobenioidCore P) {A B : C}
    (hbt : IsBaseTrivial P A) (haa : IsAutAmple P A)
    (φD : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base) :
    ∃ χ : A ⟶ B, IsPullBack P χ ∧ P.Base χ = φD := by
  obtain ⟨X, ψ, hψpb, θ, hθ⟩ := plBk_realize P F B φD
  exact prop_1_11_i_lift P F hbt haa φD ψ hψpb θ hθ

include P in
/-- ★★**(vi) の mono の `⟹` 向き** —— pull-back 射が mono なら底も mono。

★★**`plBk_realize`(圏同値)と `Definition 1.2, (ii)` の両方が要る**:
1. `u : Y ⟶ A_𝒟` を pull-back 射 `ψ : X ⟶ A` の底として実現する(`plBk_realize`)
2. `v` については、★**`Definition 1.2, (ii)` の全射性**で `f : X ⟶ A`
   (底が `θ ≫ v`、`f ≫ φ = ψ ≫ φ`)を作る
   ——★**「`ψ ≫ φ` を第1成分に取る」のが要点**である
3. `Mono φ` から `f = ψ`、よって底が一致し、`θ` が同型なので `u = v`

★**原文が「together with」で並べた 2 つの道具が、ここでは 2 と 1 に対応する。** -/
theorem prop_1_11_vi_baseMono_of_mono (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hφ : IsPullBack P φ) (hm : Mono φ) : Mono (P.Base φ) := by
  refine ⟨fun {Y} u v h => ?_⟩
  -- 手1: `u` を pull-back 射の底として実現
  obtain ⟨X, ψ, hψpb, θ, hθ⟩ := plBk_realize P F A u
  -- 手2: `v` 側の持ち上げを `Definition 1.2, (ii)` の全射性で作る
  obtain ⟨_, hsurj⟩ := hφ X
  obtain ⟨f, hf⟩ := hsurj ⟨(ψ ≫ φ, θ.hom ≫ v), by
    rw [P.Base_comp, hθ, Category.assoc, Category.assoc, h]⟩
  have hf' := congrArg Subtype.val hf
  -- 手3: `Mono φ` から `f = ψ`
  have hfψ : f = ψ := hm.right_cancellation _ _ (congrArg Prod.fst hf')
  have hbv : P.Base f = θ.hom ≫ v := congrArg Prod.snd hf'
  have hbase : θ.hom ≫ u = θ.hom ≫ v := by
    rw [← hθ, ← hbv, hfψ]
  haveI : IsIso θ.hom := inferInstance
  exact (cancel_epi θ.hom).mp hbase

include P in
/-- ★**(vi) の mono の完成形**(iff)。 -/
theorem prop_1_11_vi_mono (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hφ : IsPullBack P φ) : Mono φ ↔ Mono (P.Base φ) :=
  ⟨prop_1_11_vi_baseMono_of_mono P F φ hφ, prop_1_11_vi_mono_of_baseMono P φ hφ⟩

include P in
/-- ★★**(vi) の fiberwise-surjective の `⟸` 向き**。

★**手順**: 底で得た `δZ_𝒟` を `plBk_realize` で `𝒞` へ持ち上げ、
★**`Definition 1.2, (ii)` の全射性で `δA` を作る**(第1成分に `ψ ≫ γ` を取る)。 -/
theorem prop_1_11_vi_fs_of_baseFs (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hφ : IsPullBack P φ) (hbfs : IsFiberwiseSurjective (P.Base φ)) :
    IsFiberwiseSurjective φ := by
  intro Z γ
  obtain ⟨DD, δAD, δZD, hsq⟩ := hbfs (P.Base γ)
  obtain ⟨X, ψ, hψpb, θ, hθ⟩ := plBk_realize P F Z δZD
  obtain ⟨_, hsurj⟩ := hφ X
  obtain ⟨δA, hδA⟩ := hsurj ⟨(ψ ≫ γ, θ.hom ≫ δAD), by
    rw [P.Base_comp, hθ, Category.assoc, Category.assoc, hsq]⟩
  exact ⟨X, δA, ψ, congrArg Prod.fst (congrArg Subtype.val hδA)⟩

include P in
/-- ★★**(vi) の fiberwise-surjective の `⟹` 向き**。

★**手順**: 底の射 `γ_𝒟` を `plBk_realize` で `𝒞` へ持ち上げ、
`𝒞` 側の fiberwise-surjectivity を当てて、底を取る。 -/
theorem prop_1_11_vi_baseFs_of_fs (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hfs : IsFiberwiseSurjective φ) :
    IsFiberwiseSurjective (P.Base φ) := by
  intro ZD γD
  obtain ⟨X, ψ, hψpb, θ, hθ⟩ := plBk_realize P F B γD
  obtain ⟨Dd, δA, δX, hsq⟩ := hfs ψ
  refine ⟨(P.toElem.obj Dd).base, P.Base δA, P.Base δX ≫ θ.hom, ?_⟩
  have h := congrArg P.Base hsq
  rw [P.Base_comp, P.Base_comp, hθ] at h
  rw [h, Category.assoc]

include P in
/-- ★★**(vi) の fiberwise-surjective の完成形**(iff)。 -/
theorem prop_1_11_vi_fs (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hφ : IsPullBack P φ) :
    IsFiberwiseSurjective φ ↔ IsFiberwiseSurjective (P.Base φ) :=
  ⟨prop_1_11_vi_baseFs_of_fs P F φ, prop_1_11_vi_fs_of_baseFs P F φ hφ⟩

include P in
/-- ★★**(vi) の FSM の完成形**(iff)。

★**`IsFSMMorphism = IsFiberwiseSurjective ∧ Mono`** なので、
上の 2 つの iff を組むだけ。★**4 つのうち 3 つがこれで揃った。** -/
theorem prop_1_11_vi_fsm (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hφ : IsPullBack P φ) :
    IsFSMMorphism φ ↔ IsFSMMorphism (P.Base φ) := by
  constructor
  · rintro ⟨hfs, hm⟩
    exact ⟨(prop_1_11_vi_fs P F φ hφ).mp hfs, (prop_1_11_vi_mono P F φ hφ).mp hm⟩
  · rintro ⟨hfs, hm⟩
    exact ⟨(prop_1_11_vi_fs P F φ hφ).mpr hfs, (prop_1_11_vi_mono P F φ hφ).mpr hm⟩

include P in
/-- ★★**(vi) の irreducible の `⟹` 向き** —— pull-back 射が irreducible なら底も。

★**手順**(mono・fiberwise-surjective と同じ骨格):
1. 底の分解 `Base φ = v ≫ u` の `u` を `plBk_realize` で `𝒞` へ持ち上げ `ψ` を得る
2. ★**`Definition 1.2, (ii)` の全射性**で `β : A ⟶ X`(`β ≫ ψ = φ`)を作る
   ——★**第1成分に `φ` 自身を取る**のが今回の工夫
3. `𝒞` 側の irreducibility を使い、どちらが同型かで場合分けして底へ落とす -/
theorem prop_1_11_vi_baseIrred_of_irred (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hφ : IsPullBack P φ) (hirr : IsIrreducibleMor φ) :
    IsIrreducibleMor (P.Base φ) := by
  constructor
  · intro hiso
    exact hirr.1 (isIso_of_isPullBack_of_isBaseIso P F φ hφ hiso)
  · intro Y v u hfac
    obtain ⟨X, ψ, hψpb, θ, hθ⟩ := plBk_realize P F B u
    haveI : IsIso θ.hom := inferInstance
    haveI : IsIso θ.inv := inferInstance
    obtain ⟨_, hsurj⟩ := hψpb A
    obtain ⟨β, hβ⟩ := hsurj ⟨(φ, v ≫ θ.inv), by
      rw [hfac, hθ, Category.assoc, ← Category.assoc θ.inv, θ.inv_hom_id,
        Category.id_comp]⟩
    have hβ' := congrArg Subtype.val hβ
    have hβψ : β ≫ ψ = φ := congrArg Prod.fst hβ'
    have hβb : P.Base β = v ≫ θ.inv := congrArg Prod.snd hβ'
    rcases hirr.2 X β ψ hβψ.symm with h | h
    · left
      haveI := h
      haveI : IsIso (P.Base β) := by
        refine ⟨P.Base (inv β), ?_, ?_⟩
        · rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id]
        · rw [← P.Base_comp, IsIso.inv_hom_id, P.Base_id]
      rw [hβb] at this
      exact IsIso.of_isIso_comp_right v θ.inv
    · right
      haveI := h
      haveI : IsIso (P.Base ψ) := by
        refine ⟨P.Base (inv ψ), ?_, ?_⟩
        · rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id]
        · rw [← P.Base_comp, IsIso.inv_hom_id, P.Base_id]
      rw [hθ] at this
      exact IsIso.of_isIso_comp_left θ.hom u

include P in
/-- ★★**(vi) の irreducible の `⟸` 向き** —— 底が irreducible なら pull-back 射も。

★★**要点は `Proposition 1.7, (v)`**: **pull-back 射の両因子は pull-back 射である。**
だから底で同型が出た側を `isIso_of_isPullBack_of_isBaseIso` で `𝒞` へ上げられる。

★**`⟹` 向きと違い、圏同値は要らない** —— `Proposition 1.7, (v)` と
`Proposition 1.9` の「底が同型な pull-back は同型」だけで出る。
★**同じ iff の 2 つの向きが、まったく違う道具を使う。** -/
theorem prop_1_11_vi_irred_of_baseIrred (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hφ : IsPullBack P φ) (hbirr : IsIrreducibleMor (P.Base φ)) :
    IsIrreducibleMor φ := by
  constructor
  · intro hiso
    haveI := hiso
    refine hbirr.1 ⟨P.Base (inv φ), ?_, ?_⟩
    · rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id]
    · rw [← P.Base_comp, IsIso.inv_hom_id, P.Base_id]
  · intro X β α hfac
    -- ★両因子が pull-back(`Proposition 1.7, (v)`)
    obtain ⟨hβpb, hαpb⟩ := prop_1_7_v_pullBack P F β α (hfac ▸ hφ)
    have hb : P.Base φ = P.Base β ≫ P.Base α := by rw [hfac, P.Base_comp]
    rcases hbirr.2 _ (P.Base β) (P.Base α) hb with h | h
    · exact Or.inl (isIso_of_isPullBack_of_isBaseIso P F β hβpb h)
    · exact Or.inr (isIso_of_isPullBack_of_isBaseIso P F α hαpb h)

include P in
/-- ★★★**(vi) の irreducible の完成形**(iff)。★**これで (vi) の 4 性質がすべて揃った。** -/
theorem prop_1_11_vi_irred (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hφ : IsPullBack P φ) :
    IsIrreducibleMor φ ↔ IsIrreducibleMor (P.Base φ) :=
  ⟨prop_1_11_vi_baseIrred_of_irred P F φ hφ, prop_1_11_vi_irred_of_baseIrred P F φ hφ⟩

/-! ## ★(vii) の本体 —— 「WLOG」を合成閉性として明示する

★**原文**(p.38、目視):
> By applying the factorizations of Definition 1.3, (iv), (a); Definition 1.3, (v), (b),
> it follows **immediately** that we may assume **without loss of generality** …
> that ϵ is a pull-back morphism, an isometric pre-step, a co-angular pre-step,
> or a morphism of Frobenius type.

★★**「without loss of generality」の中身は合成閉性である。**
任意の `ϵ` は 4 種類の合成に分解される(`Definition 1.3, (iv), (a)` と `(v), (b)`)ので、
★**「各種類で成り立つ」＋「合成で閉じる」⟹「任意の `ϵ` で成り立つ」。**

★**原文は「immediately」と書くが、合成閉性は 2 回の適用を要し、
しかも 2 回目は 1 回目が返した `γ` を入力にする。**
★**`Proposition 1.10, (i)` の「バケツリレー」と同じ形**である。
-/

/-- ★(vii) が主張する性質。★**「WLOG」を扱うには名前が要る。**

`ε : Cc ⟶ B` について、**`B` へ入る任意の co-angular pre-step `φ` を、
`ε` に沿って引き戻せる**こと。 -/
def LiftsCoaPre {B Cc : C} (ε : Cc ⟶ B) : Prop :=
  ∀ {A : C} (φ : A ⟶ B), IsCoAngular P φ → IsPreStep P φ →
    ∃ (Dd : C) (γ : Dd ⟶ Cc) (α : Dd ⟶ A),
      IsCoAngular P γ ∧ IsPreStep P γ ∧ γ ≫ ε = α ≫ φ

include P in
/-- ★★**(vii) の「WLOG」の実体** —— 性質は合成で閉じる。

★**2 回目の適用が 1 回目の出力 `γ₂` を入力にする**のが要点。
★**原文の「immediately」は、この受け渡しを含んでいる。** -/
theorem liftsCoaPre_comp {B E Cc : C} {ε₁ : Cc ⟶ E} {ε₂ : E ⟶ B}
    (h₁ : LiftsCoaPre P ε₁) (h₂ : LiftsCoaPre P ε₂) : LiftsCoaPre P (ε₁ ≫ ε₂) := by
  intro A φ hφc hφs
  obtain ⟨D₂, γ₂, α₂, hγ₂c, hγ₂s, hsq₂⟩ := h₂ φ hφc hφs
  obtain ⟨D₁, γ₁, α₁, hγ₁c, hγ₁s, hsq₁⟩ := h₁ γ₂ hγ₂c hγ₂s
  refine ⟨D₁, γ₁, α₁ ≫ α₂, hγ₁c, hγ₁s, ?_⟩
  rw [← Category.assoc, hsq₁, Category.assoc, hsq₂, ← Category.assoc]

include P in
/-- ★**同型については自明に成り立つ**(基底の場合)。

★`γ := φ ≫ inv ε`、`α := 𝟙` に取ればよい。 -/
theorem liftsCoaPre_of_isIso (F : FrobenioidCore P) {B Cc : C} (ε : Cc ⟶ B) [IsIso ε] :
    LiftsCoaPre P ε := by
  intro A φ hφc hφs
  refine ⟨A, φ ≫ inv ε, 𝟙 A, ?_, ?_, ?_⟩
  · exact F.coAngularComp φ (inv ε) hφc (isCoAngular_of_isIso P (inv ε))
  · exact IsPreStep.comp P hφs (isPreStep_of_isIso P (inv ε))
  · rw [Category.assoc, IsIso.inv_hom_id, Category.comp_id, Category.id_comp]

/-! ### ★スライス側の実現補題 —— `coaPre_realize` の対応物

★`Proposition 1.10` で作った `coaPre_realize` は**コスライス**(対象から出る射)側だった。
(vii) の co-angular pre-step の場合には**スライス**(対象へ入る射)側が要る。

★★**`Order` は前順序圏なので、`op` の同型は両向きの `MLe`** であり、
`mle_antisymm`(integral ＋ sharp)で等号になる ——★**`coaPre_realize` と同じ構造。**
-/

include P in
/-- ★★**第2の圏同値のスライス側の消費** —— `Φ(B_𝒟)` の任意の元は、
**`B` へ入る co-angular pre-step の不変量**として実現できる。 -/
theorem coaPre_realize_over (G : Frobenioid P) (B : C)
    (c : Φ.val (P.toElem.obj B).base) :
    ∃ (Dd : C) (δ : Dd ⟶ B) (hδc : IsCoAngular P δ) (hδs : IsPreStep P δ),
      haveI : IsIso (P.Base δ) := hδs.2
      Φ.map (inv (P.Base δ)) (P.Div δ) = c := by
  letI := coaPreProp_isMultiplicative P G.core.coAngularComp
  haveI := G.coaPreOverEquiv B
  obtain ⟨Z, ⟨e⟩⟩ := Functor.EssSurj.mem_essImage (F := coaPreOverFunctor P B)
    (Opposite.op (toOrderCat c))
  refine ⟨Z.left.obj, Z.hom.hom, Z.hom.property.1, Z.hom.property.2, ?_⟩
  refine mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2 ?_ ?_
  · exact leOfHom e.inv.unop
  · exact leOfHom e.hom.unop

include P in
/-- ★★★**(vii) の co-angular pre-step の場合** ——
原文の「follows immediately from the second equivalence of categories of
Definition 1.3, (iii), (d)」の実体。

★★**構成**: `ε` と `φ` の不変量 `x_ε`、`x_φ` の**和**を取り、
それを不変量とする co-angular pre-step `δ : Dd ⟶ B` を実現する。
`x_ε ≤ x_ε + x_φ` と `x_φ ≤ x_ε + x_φ` から、橋(`coaPre_factor_of_mle`)が
`γ : Dd ⟶ Cc` と `α : Dd ⟶ A` を与え、どちらも `δ` を経由するので等しくなる。

★★**「和を取る」が要点である。** 原文は「immediately」と書くが、
★**`Order(Φ(B))` の中で `x_ε` と `x_φ` の上界を作る**という一手が入っている。
★**モノイドが可換だから和が上界になる** —— そこも書かれていない。 -/
theorem prop_1_11_vii_coaPre (G : Frobenioid P) {B Cc : C}
    (ε : Cc ⟶ B) (hεc : IsCoAngular P ε) (hεs : IsPreStep P ε) :
    LiftsCoaPre P ε := by
  intro A φ hφc hφs
  haveI hbε : IsIso (P.Base ε) := hεs.2
  haveI hbφ : IsIso (P.Base φ) := hφs.2
  set xε := Φ.map (inv (P.Base ε)) (P.Div ε) with hxε
  set xφ := Φ.map (inv (P.Base φ)) (P.Div φ) with hxφ
  obtain ⟨Dd, δ, hδc, hδs, hδinv⟩ := coaPre_realize_over P G B (xε + xφ)
  haveI hbδ : IsIso (P.Base δ) := hδs.2
  have hleε : MLe xε (Φ.map (inv (P.Base δ)) (P.Div δ)) := by
    rw [hδinv]; exact ⟨xφ, rfl⟩
  have hleφ : MLe xφ (Φ.map (inv (P.Base δ)) (P.Div δ)) := by
    rw [hδinv]; exact ⟨xε, by rw [add_comm]⟩
  obtain ⟨γ, hγc, hγs, hγ⟩ := coaPre_factor_of_mle P G ε hεc hεs δ hδc hδs hleε
  obtain ⟨α, _, _, hα⟩ := coaPre_factor_of_mle P G φ hφc hφs δ hδc hδs hleφ
  exact ⟨Dd, γ, α, hγc, hγs, by rw [hγ, hα]⟩

/-! ## ★(v) —— 圏同値の関手性

原文 (FrdI p.37):
> (v) The equivalences of categories of Definition 1.3, (iii), (d), are “functorial”

原文 (FrdI p.37):
> if and only if ψ is.

★**原文の証明**(p.38、目視):
> the uniqueness of ψ follows from the fact that **β is a monomorphism**
> [cf. Definition 1.3, (v), (a)] in the non-resp'd case and from the
> **total epimorphicity of C applied to α** in the resp'd case.

★★**2 つの向きで一意性の理由が違う** —— non-resp'd は **`β` が mono**、
resp'd は **`α` が epi**。
★**`Proposition 1.11, (iv)` でも同じことが起きた**(単射性は epi、一意性は mono)。
★★**この命題は「向き」に一貫して自覚的である。**
-/

include P in
/-- ★**(v) の一意性(non-resp'd)** —— `β` が mono であることから。 -/
theorem prop_1_11_v_unique_nonresp (F : FrobenioidCore P) {A B Cc Dd : C}
    (φ : A ⟶ B) (α : Cc ⟶ A) (β : Dd ⟶ B) (hβs : IsPreStep P β)
    {ψ₁ ψ₂ : Cc ⟶ Dd}
    (h₁ : ψ₁ ≫ β = α ≫ φ) (h₂ : ψ₂ ≫ β = α ≫ φ) : ψ₁ = ψ₂ := by
  haveI : Mono β := F.preStepMono β hβs
  exact (cancel_mono β).mp (h₁.trans h₂.symm)

include P in
/-- ★**(v) の一意性(resp'd)** —— `𝒞` の totally epimorphicity から。

★★**non-resp'd とは理由が違う** —— あちらは mono、こちらは epi。
★**原文はこの違いを明示している。** -/
theorem prop_1_11_v_unique_resp {A B Cc Dd : C}
    (φ : A ⟶ B) (α : A ⟶ Cc) (β : B ⟶ Dd)
    {ψ₁ ψ₂ : Cc ⟶ Dd}
    (h₁ : α ≫ ψ₁ = φ ≫ β) (h₂ : α ≫ ψ₂ = φ ≫ β) : ψ₁ = ψ₂ := by
  haveI : Epi α := P.totEpiC _ _ _
  exact (cancel_epi α).mp (h₁.trans h₂.symm)

include P in
/-- ★★**(v) の存在(non-resp'd、`φ` が co-angular pre-step の場合)**。

原文 (FrdI p.38):
> When φ is a co-angular pre-step, the existence of a co-angular pre-step ψ as desired
> follows formally from the equivalences of categories of Definition 1.3, (iii), (d).

★★**「formally」の中身は不変量の計算である**:
`α ≫ φ` の不変量を計算すると ★**`x_φ + x_β`** になる
(仮定 `x_α = Φ.map (Base φ) x_β` を使う)。
したがって `MLe x_β (x_φ + x_β)` が成り立ち、橋(`coaPre_factor_of_mle`)が `ψ` を与える。

★**「和が上界」がここでも効く** —— (vii) の co-angular pre-step の場合と同じ形。 -/
theorem prop_1_11_v_exists_nonresp (G : Frobenioid P) {A B Cc Dd : C}
    (φ : A ⟶ B) (hφc : IsCoAngular P φ) (hφs : IsPreStep P φ)
    (α : Cc ⟶ A) (hαc : IsCoAngular P α) (hαs : IsPreStep P α)
    (β : Dd ⟶ B) (hβc : IsCoAngular P β) (hβs : IsPreStep P β)
    (hcond : haveI : IsIso (P.Base α) := hαs.2
             haveI : IsIso (P.Base β) := hβs.2
             Φ.map (inv (P.Base α)) (P.Div α)
               = Φ.map (P.Base φ) (Φ.map (inv (P.Base β)) (P.Div β))) :
    ∃ ψ : Cc ⟶ Dd, IsCoAngular P ψ ∧ IsPreStep P ψ ∧ ψ ≫ β = α ≫ φ := by
  haveI hbα : IsIso (P.Base α) := hαs.2
  haveI hbβ : IsIso (P.Base β) := hβs.2
  haveI hbφ : IsIso (P.Base φ) := hφs.2
  have hcomp_c : IsCoAngular P (α ≫ φ) := G.core.coAngularComp α φ hαc hφc
  have hcomp_s : IsPreStep P (α ≫ φ) := IsPreStep.comp P hαs hφs
  haveI hbc : IsIso (P.Base (α ≫ φ)) := hcomp_s.2
  -- ★不変量の計算: `inv (α ≫ φ)` の不変量 = `x_φ + x_β`
  have hinv : Φ.map (inv (P.Base (α ≫ φ))) (P.Div (α ≫ φ))
      = Φ.map (inv (P.Base φ)) (P.Div φ) + Φ.map (inv (P.Base β)) (P.Div β) := by
    have hbase : P.Base (α ≫ φ) = P.Base α ≫ P.Base φ := P.Base_comp _ _
    have hdiv : P.Div (α ≫ φ)
        = Φ.map (P.Base α) (P.Div φ) + P.Div α := by
      rw [P.Div_comp, show P.degFr φ = 1 from hφs.1]
      simp
    have hi : inv (P.Base (α ≫ φ)) = inv (P.Base φ) ≫ inv (P.Base α) := by
      refine CategoryTheory.IsIso.inv_eq_of_hom_inv_id ?_
      rw [hbase, Category.assoc, ← Category.assoc (P.Base φ), IsIso.hom_inv_id,
        Category.id_comp, IsIso.hom_inv_id]
    rw [hi, hdiv, Φ.map_comp, map_add, ← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id,
      map_add, hcond, ← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id]
  refine coaPre_factor_of_mle P G β hβc hβs (α ≫ φ) hcomp_c hcomp_s ?_
  rw [hinv]
  exact ⟨Φ.map (inv (P.Base φ)) (P.Div φ), by rw [add_comm]⟩

include P in
/-- ★**(v) の「pull-back ⟺」の還元** —— 両方が linear なら、
「pull-back ⟺」は「**LB-invertible ⟺**」に落ちる。

★**根拠は `Proposition 1.4, (ii)`**: `IsPullBack φ ⟺ IsLBInvertible φ ∧ IsLinear φ`。
`𝒞^lin` の射は定義から linear なので、**linear の部分は自動**であり、
★**残るのは LB-invertible の部分だけ**である。

★★**原文は「Moreover, φ is a pull-back morphism if and only if ψ is」と
1 文で書くが、`𝒞^lin` の中で言っている以上、実際に比べているのは
LB-invertible 性だけである。** ★**この還元は原文に書かれていない。** -/
theorem prop_1_11_v_pullBack_iff_reduce (F : FrobenioidCore P) {A B Cc Dd : C}
    (φ : A ⟶ B) (ψ : Cc ⟶ Dd) (hφl : IsLinear P φ) (hψl : IsLinear P ψ) :
    (IsPullBack P φ ↔ IsPullBack P ψ) ↔ (IsLBInvertible P φ ↔ IsLBInvertible P ψ) := by
  constructor
  · intro h
    constructor
    · intro hlb
      exact ((prop_1_4_ii P F ψ).mp (h.mp ((prop_1_4_ii P F φ).mpr ⟨hlb, hφl⟩))).1
    · intro hlb
      exact ((prop_1_4_ii P F φ).mp (h.mpr ((prop_1_4_ii P F ψ).mpr ⟨hlb, hψl⟩))).1
  · intro h
    constructor
    · intro hpb
      exact (prop_1_4_ii P F ψ).mpr ⟨h.mp ((prop_1_4_ii P F φ).mp hpb).1, hψl⟩
    · intro hpb
      exact (prop_1_4_ii P F φ).mpr ⟨h.mpr ((prop_1_4_ii P F ψ).mp hpb).1, hφl⟩

/-! ### ★★残りを型で書く —— (vii) の pull-back の場合は (v) を経由する

★**原文**(p.38、目視):
> If ϵ is a pull-back morphism, then it follows immediately
> [**by "pulling back the zero divisor of φ via ϵ"** — cf. assertion (v)]
> that there exist a pull-back morphism α : D →A and a co-angular pre-step γ : D →C
> such that ϵ ◦γ = φ ◦α.

★★**「pulling back the zero divisor of φ via ϵ」の意味が特定できた**:
`φ` の不変量 `x_φ ∈ Φ(B_𝒟)` を `Φ.map (Base ε)` で `Φ(Cc_𝒟)` へ引き戻し、
★**それを不変量とする co-angular pre-step `γ : Dd ⟶ Cc` を実現する**
(`coaPre_realize_over` —— 我々が作った道具)。

★**そのうえで (v) を `φ_{(v)} := ε`、`α_{(v)} := γ`、`β_{(v)} := φ` に当てる**と
`ψ : Dd ⟶ A` が出て、それが求める `α` である。

★★**したがって (vii) の pull-back の場合は、(v) の存在を
「`φ` が pull-back の場合」について要求する。**
我々が実装した (v) の存在は「`φ` が co-angular pre-step の場合」だけなので、
★**そこが穴である。**

★**原文が (v) の pull-back の場合に挙げる道具は 4 つ**:
`Definition 1.3, (i), (c)` の圏同値 / `Definition 1.2, (ii)` /
`Proposition 1.7, (i), (v)` / `Definition 1.3, (iii), (d)` の圏同値。
★**4 つとも我々は持っている**(`plBk_realize` / `IsPullBack` /
`prop_1_7_v_pullBack` / `coaPre_factor_of_mle`)。
★**足りないのは組み方であって、材料ではない。**

★**「試していない」と明記する**(3 分類の第 3 カテゴリ)。
★**ただし穴の位置は型で確定している**: (v) の存在を `φ` が pull-back の場合に拡張すること。
-/

/-! ### ★★(v) の pull-back の場合を試した —— なぜ co-angular の手が移らないか

前段で「材料は全部あるが組み方が分からない」と書いた。★**試した。** 結果:

★★**co-angular pre-step の場合の手は、そのままでは移らない。理由が特定できた。**

co-angular の場合、我々は
> `α ≫ φ` の不変量を計算し、`coaPre_factor_of_mle`(スライスの圏同値)を当てる
という手を使った。★**これには `α ≫ φ` が co-angular pre-step であることが要る。**

★**`φ` が pull-back のときは、`α ≫ φ` は pre-step ではない**:
- `α` は co-angular pre-step(linear ＋ base-isomorphism)
- `φ` は pull-back(co-angular ＋ linear だが、★**base-isomorphism とは限らない**)
- よって `α ≫ φ` は **co-angular ＋ linear** だが **base-isomorphism ではない**

★★**`(𝒞^coa-pre)_B` の対象になれないので、スライスの圏同値が使えない。**

★**これが「co-angular の場合」と「pull-back の場合」で
原文が別々の道具を挙げている理由**である:
- co-angular pre-step の場合 → `Definition 1.3, (iii), (d)` の圏同値**だけ**
- pull-back の場合 → ★**4 つの道具**(`Def 1.3, (i), (c)` の圏同値 / `Def 1.2, (ii)` /
  `Prop 1.7, (i), (v)` / `Def 1.3, (iii), (d)` の圏同値)

★**原文は 4 つを並べるだけだが、なぜ 4 つ要るかは
「`α ≫ φ` が pre-step でなくなる」という一点に帰着する。**

★★**「試していない」から「試して、なぜ難しいかが分かった」へ移った**
(3 分類の第 3 → 第 2 カテゴリ)。★**この移動は 3 回目である**
(`Proposition 1.10, (iii)` の moreover で 2 回、ここで 1 回)。
-/

include P in
/-- ★★★**(v) の存在(non-resp'd、`φ` が pull-back の場合)** —— 前段で「組み方が分からない」
と書いたものが組めた。

★★**構成(5 手)**:
1. ★**`plBk_realize`** で `Base φ ≫ inv (Base β)` を底とする pull-back `φ' : E ⟶ Dd` を得る
2. ★**`Definition 1.2, (ii)` の全射性**で `β' : E ⟶ A`(`β' ≫ φ = φ' ≫ β`、底は `θ`)を作る
3. `β'` が co-angular pre-step であることを **`Proposition 1.7, (v)`** から出す
4. ★**`β'` の不変量が `Φ.map (Base φ) x_β` に等しい**ことを計算する
   ——★**仮定より `x_α` に等しい**
5. 不変量が等しいので橋(`coaPre_factor_of_mle`)が `κ : Cc ⟶ E`(`κ ≫ β' = α`)を与え、
   ★**`ψ := κ ≫ φ'`**

★★**原文が挙げる 4 つの道具が、手 1・2・3・5 に正確に対応する。**
★**「材料は全部ある。足りないのは組み方」と前段で書いたが、
組み方は「不変量が一致する対象を作る」だった。** -/
theorem prop_1_11_v_exists_pullBack (G : Frobenioid P) {A B Cc Dd : C}
    (φ : A ⟶ B) (hφpb : IsPullBack P φ)
    (α : Cc ⟶ A) (hαc : IsCoAngular P α) (hαs : IsPreStep P α)
    (β : Dd ⟶ B) (hβc : IsCoAngular P β) (hβs : IsPreStep P β)
    (hcond : haveI : IsIso (P.Base α) := hαs.2
             haveI : IsIso (P.Base β) := hβs.2
             Φ.map (inv (P.Base α)) (P.Div α)
               = Φ.map (P.Base φ) (Φ.map (inv (P.Base β)) (P.Div β))) :
    ∃ ψ : Cc ⟶ Dd, ψ ≫ β = α ≫ φ := by
  haveI hbα : IsIso (P.Base α) := hαs.2
  haveI hbβ : IsIso (P.Base β) := hβs.2
  obtain ⟨hφlb, hφlin⟩ := (prop_1_4_ii P G.core φ).mp hφpb
  -- 手1
  obtain ⟨E, φ', hφ'pb, θ, hθ⟩ := plBk_realize P G.core Dd (P.Base φ ≫ inv (P.Base β))
  obtain ⟨hφ'lb, hφ'lin⟩ := (prop_1_4_ii P G.core φ').mp hφ'pb
  -- 手2
  obtain ⟨_, hsurj⟩ := hφpb E
  obtain ⟨β', hβ'⟩ := hsurj ⟨(φ' ≫ β, θ.hom), by
    rw [P.Base_comp, hθ, Category.assoc, Category.assoc, IsIso.inv_hom_id,
      Category.comp_id]⟩
  have hβ'' := congrArg Subtype.val hβ'
  have hsq : β' ≫ φ = φ' ≫ β := congrArg Prod.fst hβ''
  have hbb : P.Base β' = θ.hom := congrArg Prod.snd hβ''
  -- 手3: `β'` は co-angular pre-step
  have hcocomp : IsCoAngular P (β' ≫ φ) := by
    rw [hsq]; exact G.core.coAngularComp φ' β hφ'lb.1 hβc
  have hlincomp : IsLinear P (β' ≫ φ) := by
    rw [hsq]; exact IsLinear.comp P hφ'lin hβs.1
  obtain ⟨⟨hβ'co, hβ'lin⟩, _⟩ :=
    prop_1_7_v_coAngularLinear P G.core β' φ hcocomp hlincomp
  haveI hbβ' : IsIso (P.Base β') := by rw [hbb]; infer_instance
  have hβ's : IsPreStep P β' := ⟨hβ'lin, hbβ'⟩
  -- 手4: 不変量が一致する
  have hdiv' : P.Div β' = Φ.map (P.Base φ') (P.Div β) := by
    have h := congrArg P.Div hsq
    rw [P.Div_comp, P.Div_comp, show P.Div φ = 0 from hφlb.2,
      show P.Div φ' = 0 from hφ'lb.2, show P.degFr φ = 1 from hφlin] at h
    simpa using h
  -- ★`inv` の下を書き換えると motive エラー(表 #2)。外から特徴づける
  have hinvb' : inv (P.Base β') = θ.inv :=
    CategoryTheory.IsIso.inv_eq_of_hom_inv_id (by rw [hbb, θ.hom_inv_id])
  have hinvb : Φ.map (inv (P.Base β')) (P.Div β')
      = Φ.map (inv (P.Base α)) (P.Div α) := by
    rw [hdiv', hcond, hinvb', ← Φ.map_comp, ← Φ.map_comp, hθ,
      ← Category.assoc, θ.inv_hom_id, Category.id_comp]
  -- 手5
  obtain ⟨κ, _, _, hκ⟩ :=
    coaPre_factor_of_mle P G β' hβ'co hβ's α hαc hαs (hinvb ▸ ⟨0, add_zero _⟩)
  exact ⟨κ ≫ φ', by rw [Category.assoc, ← hsq, ← Category.assoc, hκ]⟩

include P in
/-- ★★★**(vii) の pull-back の場合** —— 原文の
「by "pulling back the zero divisor of φ via ϵ" — cf. assertion (v)」の実体。

★★**2 手**:
1. ★**`φ` の不変量を `Φ.map (Base ε)` で引き戻し、それを不変量とする
   co-angular pre-step `γ : Dd ⟶ Cc` を `coaPre_realize_over` で実現する**
   ——★**これが「pulling back the zero divisor of φ via ϵ」である**
2. ★**(v) の pull-back の場合**を `φ_{(v)} := ε`、`α_{(v)} := γ`、`β_{(v)} := φ` に当てる

★**原文の引用符付きの語が、そのまま我々の補題の名前になった。** -/
theorem prop_1_11_vii_pullBack (G : Frobenioid P) {B Cc : C}
    (ε : Cc ⟶ B) (hεpb : IsPullBack P ε) : LiftsCoaPre P ε := by
  intro A φ hφc hφs
  haveI hbφ : IsIso (P.Base φ) := hφs.2
  -- 手1: `φ` の zero divisor を `ε` で引き戻す
  obtain ⟨Dd, γ, hγc, hγs, hγinv⟩ :=
    coaPre_realize_over P G Cc (Φ.map (P.Base ε) (Φ.map (inv (P.Base φ)) (P.Div φ)))
  -- 手2: (v) の pull-back の場合
  obtain ⟨ψ, hψ⟩ :=
    prop_1_11_v_exists_pullBack P G ε hεpb γ hγc hγs φ hφc hφs hγinv
  exact ⟨Dd, γ, ψ, hγc, hγs, hψ.symm⟩

/-! ### ★★(vii) の isometric pre-step の場合 —— なぜ別の道具が要るか

★**原文**(p.38、目視):
> if ϵ is an isometric pre-step, then the existence of γ, α with the desired properties
> follows formally from **the equivalence of categories of Proposition 1.9, (ii)** [induced by φ].

★**なぜ pull-back の場合と同じ手が使えないか、構造的な理由が特定できた。**

★★**`Proposition 1.4, (iii)`: LB-invertible な pre-step は同型である。**
`LB-invertible = co-angular ∧ isometric` なので、
★★**co-angular な isometric pre-step は同型**である。

★**したがって、同型でない isometric pre-step は決して co-angular ではない。**
★**我々の (v) の存在(co-angular pre-step の場合・pull-back の場合)はどちらも
`ε` の co-angular 性を使うので、この場合には使えない。**

★★**これが原文が isometric pre-step の場合にだけ `Proposition 1.9, (ii)` を
挙げている理由である。** ★**原文は理由を書いていないが、
「co-angular でないクラス」だから別の道具が要る、という一点に帰着する。**

★**同型の場合は既に済んでいる**(`liftsCoaPre_of_isIso`)ので、
★**残るのは「同型でない isometric pre-step」だけ**である。

## ★★使う道具を特定した（2026-08-16）

★★**`pushFunctor_isEquivalence` の本質的全射性がその道具である。**

`φ : A ⟶ B` が co-angular pre-step なら、`Proposition 1.9, (ii)` により
`pushFunctor φ` は圏同値である。
★**`ε : C ⟶ B` が isometric pre-step なら、それは
`Over (⟨B⟩ : ImtrPre P)` の対象そのもの**である。
本質的全射性から `W : Over (⟨A⟩ : ImtrPre P)` と同型
`(pushFunctor φ).obj W ≅ ⟨ε⟩` が取れ、
`pushFac : W.hom.hom ≫ φ = pushMid ≫ pushHom` が求める四角形を与える。

★★**しかし 1 点足りない**: `pushMid_spec` が与えるのは
**co-angular かつ base-isomorphism** であって、**linear ではない**。
原文の (vii) は `γ` が **co-angular pre-step**
(= linear ∧ base-isomorphism ∧ co-angular) であることを要求するので、
★**`pushMid` の linear 性を別途出す必要がある**。

★**次に当たるときの手**: `pushFac` の両辺に `degFr` を適用する。
左辺は `degFr φ * degFr (W.hom.hom) = 1`(どちらも pre-step なので linear)、
右辺は `degFr (pushHom) * degFr (pushMid)`。
★`ℕ+` の消去で両方が 1 になるはずである。
★★**つまり linear 性は `pushFac` から自動で出ると見ている。まだ試していない。**
-/

include P in
/-- ★★**co-angular な isometric pre-step は同型**。

★`Proposition 1.4, (iii)` の直接の帰結だが、★**(vii) の場合分けの構造を説明する。**
「同型でない isometric pre-step は co-angular でない」——だから
原文は isometric pre-step の場合にだけ別の道具(`Proposition 1.9, (ii)`)を挙げる。 -/
theorem isIso_of_isometric_coAngular_preStep (F : FrobenioidCore P) {A B : C}
    (ε : A ⟶ B) (hi : IsIsometric P ε) (hc : IsCoAngular P ε) (hs : IsPreStep P ε) :
    IsIso ε :=
  prop_1_4_iii P F ε ⟨hc, hi⟩ hs

/-! ### ★★(vii) の Frobenius 型の場合 —— **2 種類目の「we may assume」**

★**原文**(p.38、目視):
> Finally, we consider the case where ϵ is a morphism of Frobenius type.
> By applying the second equivalence of categories of Definition 1.3, (iii), (d),
> it follows that **we may assume [by replacing φ by the composite of φ with an
> appropriate pre-step A′ →A]** that Div(φ) = degFr(ϵ) · x, for some x ∈Φ(A).
> Thus, … there exist a morphism of Frobenius type α : D →A and a co-angular
> pre-step γ : D →C such that ϵ ◦γ = φ ◦α [**cf. also Proposition 1.10, (i)**].

★★**この証明には「we may assume」が 2 回出る。しかも種類が違う。**

| どこ | 何を仮定してよいとするか | 正当化 |
|---|---|---|
| 場合分けの前 | `ϵ` が 4 種類のどれか | ★**合成閉性**(`liftsCoaPre_comp`) |
| Frobenius 型の中 | `Div(φ) = degFr(ϵ) · x` | ★**`φ` を `ρ ≫ φ` に置き換える** |

★★**1 つ目は `ϵ` を分解する。2 つ目は `φ` を置き換える。**
★**同じ語(「we may assume」)が、反対側の射に作用している。**

★**2 つ目の正当化**: `γ ≫ ε = α′ ≫ (ρ ≫ φ)` が得られれば
`α := α′ ≫ ρ` として `γ ≫ ε = α ≫ φ` になる。★**結合則だけ。**

★**そして原文は `Proposition 1.10, (i)` を引く** ——
★★**我々がこのセッションで完成させた項目**である。
★**1.10 の投資が 1.11 の最後の場合でも回収される。**
-/

include P in
/-- ★★**(vii) の 2 種類目の「we may assume」の正当化** ——
`φ` を `ρ ≫ φ` に置き換えて解けたなら、元の `φ` でも解ける。

★**結合則だけ**だが、★**「we may assume」を型にすると
何が保証されているかが明示される。** -/
theorem liftsCoaPre_precomp {B Cc A A' : C} (ε : Cc ⟶ B) (φ : A ⟶ B) (ρ : A' ⟶ A)
    {Dd : C} (γ : Dd ⟶ Cc) (α' : Dd ⟶ A')
    (h : γ ≫ ε = α' ≫ (ρ ≫ φ)) :
    (α' ≫ ρ) ≫ φ = γ ≫ ε := by
  rw [Category.assoc, ← h]

/-! ### ★★(vii) の Frobenius 型の場合 —— 向きの制約が見えた

★**原文**(p.38、目視)は、この場合の結論を
> there exist **a morphism of Frobenius type α : D →A** and a co-angular pre-step γ : D →C
と書く。★**`α` は `A` へ**入る**Frobenius 型射である。**

★★**ところが `Definition 1.3, (ii)`(`frobDegSurj`)は
「対象から**出る**Frobenius 型射」しか与えない**:
```
frobDegSurj : ∀ (A : C) (n : ℕ+), ∃ (B : C) (φ : A ⟶ B), IsFrobeniusType φ ∧ degFr φ = n
```
★**`A` へ入る Frobenius 型射の存在は、`Definition 1.3` の公理には無い。**

★**`Definition 1.2, (iv)` の `perfect object` の条件 (a) はそれを与える**:
> every B ∈Ob(C) base-isomorphic to C appears as **the codomain of** a morphism of
> Frobenius type of Frobenius degree n

★★**しかし (vii) は perfect 型を仮定していない。**

★**したがって `α` の出どころは、原文が「by applying again the second equivalence
of categories of Definition 1.3, (iii), (d)」と書いている部分にあるはず**だが、
第2の圏同値は **co-angular pre-step** を与えるものであって
Frobenius 型射を与えるものではない。

★★**ここは「試して、道具の向きが合わないことが分かった」段階である**
(3 分類の第 2 カテゴリ)。★**穴の位置は「`A` へ入る Frobenius 型射をどこから得るか」
の 1 点に確定した。**

★**これは `Proposition 1.10, (iii)` の moreover で我々が経験した
「探していた向きが逆だった」と同じ形の可能性がある** ——
★**あのときも、必要な道具は別の向きで手元にあった。**

## ★★さらに気づいたこと: 原文の証明は主張より強い

★**(vii) の主張は `α` について何も要求していない**:
> there exists a co-angular pre-step γ : D →C and **a morphism α : D →A** such that
> ϵ ◦γ = φ ◦α.

★**ところが証明の Frobenius 型の場合では
「a morphism of **Frobenius type** α : D →A」を作っている。**

★★**原文の証明は、主張が求めるより強いものを作っている。**
★**これは無駄ではない** —— 場合分けの他の場合でも `α` の型が明示されている
(pull-back の場合は「a pull-back morphism α」)。
★★**著者は各場合で `α` の型を追跡しており、主張ではそれを捨てている。**

★**測定**: 我々は「原文が書かない一歩」を数えてきたが、
★★**ここは逆に「原文が書いているが主張には要らない情報」である。**
★**証明は主張より多くを語ることがある** ——
**そしてその余剰は、後の命題で使われるかもしれない**(未確認)。
-/

include P in
/-- ★★**(vii) の Frobenius 型の場合: `α` の次数は強制される**。

★四角形 `γ ≫ ε = α ≫ φ` の次数を取ると:
- `degFr (γ ≫ ε) = degFr ε * degFr γ = degFr ε`(`γ` は pre-step)
- `degFr (α ≫ φ) = degFr φ * degFr α = degFr α`(`φ` は pre-step)

★★**したがって `degFr α = degFr ε`。**

★**これが「原文がなぜ Frobenius 型の `α` を作るか」の説明である**:
(vii) の主張は `α` の型を要求しないが、★**次数は強制される**ので、
**次数 `degFr ε` の射を作るしかなく、その自然な供給源が `Definition 1.3, (ii)` の
Frobenius 型射**である。

★★**「主張が要求しない」ことと「自由に選べる」ことは違う。**
★**余剰に見えた情報(`α` が Frobenius 型)の半分(次数)は、実は強制されていた。** -/
theorem prop_1_11_vii_degFr_forced {A B Cc Dd : C}
    (ε : Cc ⟶ B) (φ : A ⟶ B) (γ : Dd ⟶ Cc) (α : Dd ⟶ A)
    (hγs : IsPreStep P γ) (hφs : IsPreStep P φ)
    (hsq : γ ≫ ε = α ≫ φ) : P.degFr α = P.degFr ε := by
  have h := congrArg P.degFr hsq
  rw [P.degFr_comp, P.degFr_comp, show P.degFr γ = 1 from hγs.1,
    show P.degFr φ = 1 from hφs.1] at h
  simpa using h.symm

end ABC3.Found.FrdI
