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

end ABC3.Found.FrdI
