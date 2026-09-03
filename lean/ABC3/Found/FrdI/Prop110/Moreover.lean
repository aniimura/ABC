/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.GroupTheory.MonoidLocalization.GrothendieckGroup
import Mathlib.NumberTheory.PrimeCounting
import ABC3.Found.FrdI.Prop19
import ABC3.Found.FrdI.MonoidPrime
import ABC3.Found.FrdI.Prop110.Sharp

/-!
# Prop110 —— (iii) の moreover・穴を浅くする段

☆もとの 1 枚を**ファイル内の見出し**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(iii) の「moreover」—— `𝒪^▷(A)` と `𝒪^×(A)` が perfect

原文 (FrdI p.35):
> immediately [cf. Remark 1.1.1] from the fact that A is perfect [cf. Definition 1.2,

★**主張そのものの行(「then the monoids O▷(A) and O×(A) are perfect.」)は
引用できない**(事故 #2 の再来。`▷` が layout モードで消え 15/38 文字で停止)。
★**`▷` は事故 #2 で既に記録した文字である** —— 同じ文字が別の場所で再び落ちた。
証明側の隣接行を引く。

★**まず語彙の問題**: `𝒪^▷(A)` は `Submonoid (End A)` であり**乗法**モノイドである。
我々の `IsPerfectMonoid` は**加法**モノイド用(`∀ n, Bijective (n • ·)`)。
★**乗法版が要る。** 原文は両方を同じ「perfect」で呼ぶが、
★**`Φ(A)` は加法、`𝒪^▷(A)` は乗法**で、**同じ語が 2 つの構造を指している。**

★**原文の証明**(p.35、目視):
> the fact that the monoids O▷(A) and O×(A) are perfect follows immediately from the
> fact that A is perfect [cf. Definition 1.2, (iv), **applied to the base-identity
> endomorphisms of Frobenius type of the Frobenius-trivial object A**] and Frobenius-normalized.
-/

/-- ★**乗法モノイドの perfect** —— `n` 乗が全単射。

★`IsPerfectMonoid`(加法)の乗法版。★**原文が同じ語で呼ぶ 2 つを、我々は別の述語で書く。** -/
def IsPerfectMulMonoid (M : Type*) [Monoid M] : Prop :=
  ∀ n : ℕ+, Function.Bijective (fun a : M => a ^ (n : ℕ))

/-- ★**非退化(下)**: 自明モノイドは perfect。 -/
theorem isPerfectMulMonoid_punit : IsPerfectMulMonoid PUnit :=
  fun _ => ⟨fun _ _ _ => rfl, fun _ => ⟨PUnit.unit, rfl⟩⟩

/-- ★**非退化(上)**: `Multiplicative ℕ` は perfect でない —— 2 乗写像が全射でない
(`ofAdd 1` は平方の像に入らない)。 -/
theorem not_isPerfectMulMonoid_nat : ¬ IsPerfectMulMonoid (Multiplicative ℕ) := by
  intro h
  obtain ⟨a, ha⟩ := (h 2).2 (Multiplicative.ofAdd 1)
  have h1 : Multiplicative.toAdd (a ^ (2 : ℕ)) = 1 := congrArg Multiplicative.toAdd ha
  rw [show (2 : ℕ) = 1 + 1 from rfl, pow_add, pow_one] at h1
  have h2 : Multiplicative.toAdd a + Multiplicative.toAdd a = 1 := h1
  omega

/-! ### ★★(iii) の moreover —— 残る穴を型で書く

★**構成**(原文の 3 段、目視で確認):
1. `A` が Frobenius-trivial なので、各 `n` に base-identity・Frobenius 型の
   自己射 `ζ n`(次数 `n`)がある
2. `IsPerfectObj` の (b) を `φ₁ = φ₂ = ζ n`、`ψ′ := α ∈ 𝒪^▷(A)` に当てると、
   ★**一意の pre-step `β` が `ζ n ≫ α = β ≫ ζ n` を満たす**
3. **Frobenius-normalized** を `φ := ζ n`、`α := β` に当てると
   `ζ n ≫ β^n = β ≫ ζ n = ζ n ≫ α`

★★**残る穴はただ 1 つ**: 上から `β^n = α` を結論するには
★**`ζ n` で左から消せること(mono)**が要る。
★**`Definition 1.3` は pre-step が mono であること(`preStepMono`)しか言っていない。**
`ζ n` は次数 `n` なので pre-step ではない。

★★**したがってこれが「原文の `follows immediately` が使っている、
書かれていない事実」である。** 3 通りのどれかである:
  (a) Frobenioid では Frobenius 型射は常に mono —— **原文にも我々の公理にも無い**
  (b) isotropic 型ではそうなる —— **未確認**
  (c) 別の議論で mono を回避できる —— **未発見**

★**`sorry` は置かない。** 3 通りのどれかを確かめるのが次の仕事である。
★**「試していない」と明記する**(前段で立てた 3 分類の第 3 category)。

★**語彙は先に置いた**(`IsPerfectMulMonoid`)ので、
**残りは上の 1 点だけ**である。
-/

/-! ### ★★穴が 1 段浅くなった —— mono ではなく `faithfulUpToUnits` だった

前段で親は「残る穴は `ζ n` が mono かどうか」と書いた。★**それは違った。**

★**`Definition 1.3, (vi)`(`faithfulUpToUnits`)が使える**:
`α` と `β^n` はどちらも `𝒪^▷(A)` の元(base-identity ＋ linear)であり、
- **base-equivalent**: どちらも `Base = 𝟙`
- **metrically equivalent**: `ζ n ≫ α = ζ n ≫ β^n` に `Div_comp` を当てると、
  ★**`ζ n` が Frobenius 型(`Div = 0`)かつ base-identity(`Φ.map 𝟙 = id`)**なので
  両辺が `Div α`・`Div (β^n)` に潰れる ⟹ **`Div α = Div (β^n)`**

★したがって `faithfulUpToUnits` から ★**`α = β^n ≫ u`(`u ∈ 𝒪^×(A)`)**が出る。

★★**穴は「mono か」から「単元 `u` を吸収できるか」に変わった。**
★**これは前進である**: mono は**公理に無い**が、単元の吸収は
**`𝒪^×(A)` が perfect であること**(＝ (iii) の moreover のもう半分)から出るはずで、
★**同じ条の中で閉じる可能性がある。**

★**測定**: 「残る穴」を型で書いておいたおかげで、★**穴の位置が動いたことがすぐ分かった。**
散文で「あとは mono が要る」と書いていたら、`faithfulUpToUnits` を試すきっかけが無かった。
-/

include P in
/-- ★★**`α` と `β^n` は単元の違いしかない** —— `faithfulUpToUnits` の帰結。

★前段で「mono が要る」と書いた穴が、ここまで浅くなった。 -/
theorem prop_1_10_iii_otri_upToUnit (F : FrobenioidCore P) {A : C}
    (hiso : ∀ X : C, IsIsotropic P X)
    (ζn : End A) (hζb : IsBaseIdentity P ζn) (hζf : IsFrobeniusType P ζn)
    (α β : End A) (hαb : IsBaseIdentity P α) (hαl : IsLinear P α)
    (hβb : IsBaseIdentity P β) (hβl : IsLinear P β)
    (heq : (ζn : A ⟶ A) ≫ (α : A ⟶ A) = (ζn : A ⟶ A) ≫ (β : A ⟶ A)) :
    ∃ u : End A, u ∈ OTimes P A ∧ (α : A ⟶ A) = (β : A ⟶ A) ≫ (u : A ⟶ A) := by
  -- `Div α = Div β`
  have hd : P.Div (α : A ⟶ A) = P.Div (β : A ⟶ A) := by
    have h := congrArg P.Div heq
    rw [P.Div_comp, P.Div_comp, show P.Div (ζn : A ⟶ A) = 0 from hζf.1.2,
      show P.Base (ζn : A ⟶ A) = P.Base (𝟙 A) from hζb, P.Base_id] at h
    simpa using h
  -- `α`・`β` は co-angular pre-step
  have hpre : ∀ γ : End A, IsBaseIdentity P γ → IsLinear P γ → IsPreStep P (γ : A ⟶ A) := by
    intro γ hb hl
    refine ⟨hl, ?_⟩
    show IsIso (P.Base (γ : A ⟶ A))
    rw [show P.Base (γ : A ⟶ A) = P.Base (𝟙 A) from hb, P.Base_id]
    infer_instance
  have hcoa : ∀ γ : End A, IsCoAngular P (γ : A ⟶ A) :=
    fun γ => prop_1_4_i P (γ : A ⟶ A) (fun X _ => hiso X)
  exact F.faithfulUpToUnits (α : A ⟶ A) (β : A ⟶ A)
    (show P.Base (α : A ⟶ A) = P.Base (β : A ⟶ A) by
      rw [show P.Base (α : A ⟶ A) = P.Base (𝟙 A) from hαb,
        show P.Base (β : A ⟶ A) = P.Base (𝟙 A) from hβb])
    hd (hcoa α) (hpre α hαb hαl) (hcoa β) (hpre β hβb hβl)

/-! ### ★★単元の吸収が難しい理由が構造から分かった

★**残る穴**: `α = β^n ≫ u`(`u ∈ 𝒪^×(A)`)から `α = γ^n` を出したい。
`𝒪^×(A)` も perfect なら `u = v^n` と書けるが、
★**`β^n ≫ v^n = (β ≫ v)^n` には `β` と `v` の可換性が要る。**

★★**そして `𝒪^▷(A)` は一般に可換ではない。**

原典の序文(p.1、目視):
> kind of Frobenioid "FM" is the **non-commutative** monoid given by forming the
> "semi-direct product monoid" of a given commutative monoid M with the multi-
> plicative monoid N≥1

★**`𝔽_M` は明示的に非可換**である(`M` と `ℕ≥1` の半直積)。
ただし `𝒪^▷(A)` は **linear**(次数 1)な base-identity 自己射の集まりなので、
半直積の `M` 側にあたり、**素朴には可換に見える**。

★★**しかし一般の Frobenioid ではそうとは限らない。** 我々が持っている情報は
`faithfulUpToUnits`(`Definition 1.3, (vi)`)だけで、それが言うのは
★**「base と `Div` が一致する co-angular pre-step は単元の違いしかない」** ——
すなわち ★**`𝒪^▷(A)` は単元を法としてのみ `Φ(A)` に埋まる。**

★★**したがって「単元の吸収」は、この命題の中で最も構造に近い場所にある。**
`𝒪^×(A)` が何であるかは `Definition 1.3` が直接には規定していない。

★**次の仕事(3 通り、まだ試していない)**:
  (a) `𝒪^▷(A)` の可換性を Frobenioid の公理から導く —— **未確認**
  (b) `𝒪^×(A)` が中心に入ることを示す —— **未確認**
  (c) 可換性を使わない別の議論 —— **未発見**

★**「試していない」と明記する**(3 分類の第 3 カテゴリ)。
★**ただし前段と違い、今回は「なぜ難しいか」が構造から言えている** ——
★**原典が序文で「non-commutative」と宣言している当のものに触れている。**
-/

/-! ### ★★選択肢 (a) を試した —— `𝒪^▷(A)` は**単元を法として**可換

前段で 3 通りの「まだ試していない」を挙げた。★**(a) を試した。**

★**結果**: `𝒪^▷(A)` は**単元を法として可換**である。**完全な可換性は出ない。**

★**証明**: `α, β ∈ 𝒪^▷(A)` に対し `α ≫ β` と `β ≫ α` は
- **base-equivalent**: どちらも `Base = 𝟙`
- **metrically equivalent**: `Div_comp` で
  `Div (α ≫ β) = Φ.map 𝟙 (Div β) + 1 • Div α = Div β + Div α`、
  `Div (β ≫ α) = Div α + Div β`。★**`Φ(A)` が可換なので一致する。**

したがって `faithfulUpToUnits` から ★**`α ≫ β = (β ≫ α) ≫ u`(`u ∈ 𝒪^×(A)`)**。

★★**これは「試して、部分的に出た」である**(3 分類の第 2 カテゴリへ移った)。
★**単元の吸収に必要な完全な可換性は、ここからは出ない。**

★★**測定: 可換性の破れは `Φ(A)` ではなく `𝒪^×(A)` に局在している。**
`Div` は可換に振る舞い、ずれはすべて**単元**に吸い込まれる。
★**原典が序文で「non-commutative」と言うのは半直積の `ℕ≥1` 側の話だが、
`𝒪^▷(A)` の非可換性は `𝒪^×(A)` に由来する** —— **別の出どころ**である。
-/

include P in
/-- ★★**`𝒪^▷(A)` は単元を法として可換**。

★`Div` は可換に振る舞う(`Φ(A)` が可換だから)ので、
ずれはすべて `faithfulUpToUnits` の単元に吸い込まれる。 -/
theorem otri_comm_upToUnit (F : FrobenioidCore P) {A : C}
    (hiso : ∀ X : C, IsIsotropic P X) (α β : OTri P A) :
    ∃ u : End A, u ∈ OTimes P A ∧
      ((α : End A) ≫ (β : End A)) = ((β : End A) ≫ (α : End A)) ≫ (u : A ⟶ A) := by
  have hpre : ∀ γ : OTri P A, IsPreStep P (γ : End A) := by
    intro γ
    refine ⟨γ.2.2, ?_⟩
    show IsIso (P.Base (γ : End A))
    rw [show P.Base (γ : End A) = P.Base (𝟙 A) from γ.2.1, P.Base_id]
    infer_instance
  have hcoa : ∀ f : A ⟶ A, IsCoAngular P f := fun f => prop_1_4_i P f (fun X _ => hiso X)
  have hbase : ∀ γ : OTri P A, P.Base (γ : End A) = 𝟙 _ := by
    intro γ
    rw [show P.Base (γ : End A) = P.Base (𝟙 A) from γ.2.1, P.Base_id]
  -- `Div` が一致する
  have hd : P.Div ((α : End A) ≫ (β : End A)) = P.Div ((β : End A) ≫ (α : End A)) := by
    rw [P.Div_comp, P.Div_comp, hbase α, hbase β, Φ.map_id, Φ.map_id,
      show P.degFr (β : End A) = 1 from β.2.2, show P.degFr (α : End A) = 1 from α.2.2]
    simp [add_comm]
  -- base が一致する
  have hb : P.Base ((α : End A) ≫ (β : End A)) = P.Base ((β : End A) ≫ (α : End A)) := by
    rw [P.Base_comp, P.Base_comp, hbase α, hbase β]
  exact F.faithfulUpToUnits _ _ hb hd (hcoa _)
    (IsPreStep.comp P (hpre α) (hpre β)) (hcoa _) (IsPreStep.comp P (hpre β) (hpre α))

/-! ### ★選択肢 (b) の一部 —— `β^n` が同型なら `β` も同型

★**`𝒪^×(A)` の側を先に片づけようとして出た副産物**である。

★**主張**: isotropic 型で、`β` が pre-step、`β^n` が同型なら `β` は同型。
★**証明**: `β^n` が同型 ⟹ `Div (β^n) = 0`。`β` は base-identity 型の pre-step なので
`Div (β^n) = n • Div β`(★`Div_comp` を `n` 回、`Base β = 𝟙` と `degFr β = 1` で潰す)。
★**sharp から `Div β = 0`** ⟹ `β` は isometric。isotropic から co-angular なので
LB-invertible な pre-step、`Proposition 1.4, (iii)` で**同型**。

★★**これは `𝒪^▷(A)` の中で `𝒪^×(A)` が「`n` 乗で閉じている」ことの逆向きである。**
★**単元の吸収そのものは片づかないが、「どの元が単元か」の判定は 1 つ増えた。**
-/

include P in
/-- ★`γ^k` の 3 つの不変量を一度に。★**`End A` の乗法は `x * y = y ≫ x`**
(mathlib の規約)なので、`pow_succ` を開くと `γ ≫ γ^m` になる。 -/
theorem otri_pow_invariants {A : C} (γ : End A) (hb : IsBaseIdentity P γ) (hl : IsLinear P γ) :
    ∀ k : ℕ, P.Base ((γ ^ k : End A) : A ⟶ A) = P.Base (𝟙 A) ∧
      P.degFr ((γ ^ k : End A) : A ⟶ A) = 1 ∧
      P.Div ((γ ^ k : End A) : A ⟶ A) = k • P.Div (γ : A ⟶ A) := by
  intro k
  induction k with
  | zero => exact ⟨rfl, by simpa using P.degFr_id A, by simpa using P.Div_id A⟩
  | succ m ih =>
      obtain ⟨ihb, ihd, ihv⟩ := ih
      have hmul : ((γ ^ (m + 1) : End A) : A ⟶ A) = (γ : End A) ≫ (γ ^ m : End A) := by
        rw [pow_succ]; rfl
      refine ⟨?_, ?_, ?_⟩
      · rw [hmul, P.Base_comp, ihb, P.Base_id, Category.comp_id]
        rw [show P.Base (γ : A ⟶ A) = P.Base (𝟙 A) from hb, P.Base_id]
      · rw [hmul, P.degFr_comp, ihd, hl, one_mul]
      · rw [hmul, P.Div_comp, ihv, ihd,
          show P.Base (γ : A ⟶ A) = P.Base (𝟙 A) from hb, P.Base_id, Φ.map_id]
        simp [succ_nsmul, add_comm]

include P in
/-- ★★**`β^n` が同型なら `β` も同型**(isotropic 型、`β` は `𝒪^▷(A)` の元)。

★sharp が「`n • x = 0 ⟹ x = 0`」として効く。 -/
theorem isIso_of_pow_isIso (F : FrobenioidCore P) {A : C}
    (hiso : ∀ X : C, IsIsotropic P X) (β : End A)
    (hb : IsBaseIdentity P β) (hl : IsLinear P β)
    {n : ℕ} (hn : 0 < n) (hpow : IsIso ((β ^ n : End A) : A ⟶ A)) :
    IsIso (β : A ⟶ A) := by
  haveI := hpow
  have hdiv0 : P.Div ((β ^ n : End A) : A ⟶ A) = 0 :=
    (isIsometric_of_isIso P ((β ^ n : End A) : A ⟶ A))
  have hdβ : P.Div (β : A ⟶ A) = 0 :=
    eq_zero_of_nsmul_eq_zero_of_isSharp (P.divisorial _).2 hn
      (by rw [← (otri_pow_invariants P β hb hl n).2.2]; exact hdiv0)
  have hbi : IsIso (P.Base (β : A ⟶ A)) := by
    rw [show P.Base (β : A ⟶ A) = P.Base (𝟙 A) from hb, P.Base_id]
    infer_instance
  exact prop_1_4_iii P F (β : A ⟶ A)
    ⟨prop_1_4_i P (β : A ⟶ A) (fun X _ => hiso X), hdβ⟩ ⟨hl, hbi⟩

/-! ### ★★選択肢 (c) を試した —— **公理の非対称性が穴の正体だった**

★**(c)「可換性を使わない議論」を探した。** 見つからなかったが、
★★**なぜ見つからないかが公理の形から言える。**

我々に必要なのは:
> `ζ n ≫ α = ζ n ≫ β^n`(`α, β^n ∈ 𝒪^▷(A)`)から `α = β^n` を出す

すなわち ★**「Frobenius 型射 `ζ n` に沿った `𝒪^▷` の対応が単射である」**こと。

★**`Definition 1.3` が与えるものを並べる**:

| 公理 | 沿う射 | `∃!` の位置 |
|---|---|---|
| `otriFwd` (iii)(c) | ★**co-angular pre-step** | 終域側(`β ∈ 𝒪^▷(B)` が一意) |
| `otriBwd` (iii)(c) | ★**co-angular pre-step** | 始域側(`α ∈ 𝒪^▷(A)` が一意) |
| `IsPerfectObj` (b) | ★**Frobenius 型射** | ★**始域側のみ**(`ψ` が一意) |

★★**pre-step に沿っては両側の `∃!` があり(＝全単射)、
Frobenius 型射に沿っては片側の `∃!` しかない。**

★**我々が要るのは終域側の一意性**(`ζ n ≫ α = ζ n ≫ β^n ⟹ α = β^n`)であり、
★**それは Frobenius 型射については公理に無い。**

★★**これが穴の正体である。** 「原文が書いていない一歩」でも
「原文が挙げていない依存」でもなく、★**公理の非対称性**である。

★**3 通りの結論**:
- (a) 可換性 —— 試した。**単元を法としてのみ成立**
- (b) 単元が中心 —— 試した。**一部(`β^n` 同型 ⟹ `β` 同型)のみ**
- (c) 可換性なしの議論 —— ★**試した。公理の非対称性により、この形では出ない**

★★**したがって、原文の「follows immediately」は
`Definition 1.3` に無いものを使っている**か、あるいは
★**我々が `Definition 1.3` の写しで何かを落としている**かのどちらかである。
★**後者の可能性を否定していない** —— それが次に確かめることである。
-/

/-! ### ★★★原典自身が `𝒪^×(C)` の可換性を**仮定**として置いている

★前段で「我々が `Definition 1.3` の写しで何かを落としているかもしれない」と書いた。
★**原典 p.22–23 を目視で読み直した。結果:**

**(1) 我々の写しは正しい。**
- `perfect object`: 原文 `ψ′ ◦ φ₁ = φ₂ ◦ ψ` ＝ Lean の `φ₁ ≫ ψ′ = ψ ≫ φ₂` ✓
- `Frobenius-normalized`: 原文 `α^d ◦ φ = φ ◦ α` ＝ Lean の `φ ≫ α^d = α ≫ φ` ✓
- `group-like object`: 原文 `Φ(C) = 0` ✓(我々の `IsGroupLike` はその同値形)

**(2) ★★しかし別の定義に決定的なものがあった。**

原文 (FrdI p.23):
> A Frobenius-compact object of C is defined to be an object C

原文 (FrdI p.23):
> such that O×(C) is commutative,

★**続きは引用できない**(事故 #4 の再来。同じ行の `̸=` が layout モードで消える)。
400 dpi 目視では「`O×(C)pf ≠ 0`、および `Aut_𝒞(C)` の元で `𝒪^×(C)^pf` に
`ℚ>0` の元の乗法として作用するものは実際には自明に作用する」と続く。

★★★**原典は `Frobenius-compact` の定義の中で「`𝒪^×(C)` が可換であること」を
明示的な仮定として課している。**
★**すなわち `𝒪^×(C)` は一般には可換ではない** —— 可換なら仮定に書く必要がない。

## ★★これが我々の穴の位置を確定させる

(iii) の moreover を我々は
> `α = β^n ≫ u`(`u ∈ 𝒪^×(A)`)から `α = γ^n` を出したい
まで縮めた。これには `u = v^n` と `β`・`v` の可換性が要る。

★★**原典自身が `𝒪^×` の可換性を別の場所で仮定として置いている以上、
それは `Definition 1.3` からは出ない。**

★**したがって 2 つのうちどちらかである**:
- **(α)** (iii) の moreover は、原文が (iii) で挙げていない仮定
  (`𝒪^×(A)` の可換性、あるいは Frobenius-compact 性)を使っている
- **(β)** 可換性を回避する議論があり、我々がまだ見つけていない

★**(β) を否定していない。** ただし ★**(α) の側に、原典自身が置いた証拠が 1 つある。**

★★**これは「原文が挙げる根拠が実際の根拠の一部」(4 種類目)の、
最も強い実例である** —— **原文は同じ論文の別の場所で、必要な仮定を明示している。**
-/

/-! ### ★★★選択肢 (β) が出た —— **探していた道具は最初から手元にあった**

★我々は「`ζ n` で**左から**消せること(mono)」を探し、`Definition 1.3` が
pre-step にしか mono を与えないことで詰まっていた。★**必要だったのは 2 つとも別物だった。**

★★**全射性**: `α ∈ 𝒪^▷(A)` に対し `IsPerfectObj` の (b) が pre-step `β` を与え
`ζ n ≫ α = β ≫ ζ n`。**Frobenius-normalized** から `ζ n ≫ β^n = β ≫ ζ n = ζ n ≫ α`。
★★**`𝒞` は totally epimorphic なので `ζ n` は epi** —— `Epi f` は
**`f ≫ x = f ≫ y ⟹ x = y`**、すなわち**左に置かれた `f` を消す**。
★**これがちょうど要る向きだった** ⟹ **`β^n = α`**。

★★**単射性**: `β₁^n = β₂^n = α` なら、Frobenius-normalized から
`ζ n ≫ α = βᵢ ≫ ζ n`(`i = 1, 2`)。
★★**(b) の `∃!` が「`ψ′ = α` に対する `ψ` は一意」と言っている** ⟹ **`β₁ = β₂`**。

★★★**我々は `faithfulUpToUnits` を経由して単元の吸収に苦しんでいたが、
`totEpiC` と (b) の `∃!` だけで足りた。**
★**「読み落とし」(6 種類目のギャップ)の 3 例目である。**
★★**3 回とも、必要なものは既に手元にあった。**
(1: §0 の torsion-free、2: `Gp M` の群構造の import、3: これ)

★**そして `𝒪^×(A)` の可換性は要らなかった** ——
原典が `Frobenius-compact` で可換性を仮定しているのは**別の目的**である。
★**我々が「穴」と呼んでいたものは、遠回りの経路にしか存在しなかった。**
-/

include P in
/-- ★★★**`Proposition 1.10, (iii)` の moreover(`𝒪^▷(A)` の側)完成** ——
`𝒪^▷(A)` では `n` 乗根が**一意に存在する**。

★**全射性**は `𝒞` の totally epimorphic 性(epi)から、
★**単射性**は `IsPerfectObj` の (b) の `∃!` から。
★**`𝒪^×(A)` の可換性は要らない。** -/
theorem prop_1_10_iii_otri_perfect {A : C}
    (hperf : IsPerfectObj P A) (hfn : IsFrobeniusNormalized P A)
    (n : ℕ+) (ζn : End A) (hζb : IsBaseIdentity P ζn) (hζf : IsFrobeniusType P ζn)
    (hζd : P.degFr ζn = n) (α : End A) (hα : α ∈ OTri P A) :
    ∃! β : End A, β ∈ OTri P A ∧ (β ^ (n : ℕ) : End A) = α := by
  haveI : Epi (ζn : A ⟶ A) := P.totEpiC _ _ _
  have hαs : IsPreStep P (α : A ⟶ A) := by
    refine ⟨hα.2, ?_⟩
    show IsIso (P.Base (α : A ⟶ A))
    rw [show P.Base (α : A ⟶ A) = P.Base (𝟙 A) from hα.1, P.Base_id]
    infer_instance
  have hAA : BaseIsomorphic P A A := ⟨Iso.refl _⟩
  obtain ⟨β, ⟨hβs, hsq⟩, huniq⟩ :=
    (hperf n).2 A A A A ζn ζn hζf hζd hζf hζd hAA hAA (α : End A) hαs
  -- `β` は base-identity
  have hβb : IsBaseIdentity P β := by
    show P.Base β = P.Base (𝟙 A)
    have h := congrArg P.Base hsq
    rw [P.Base_comp, P.Base_comp, show P.Base (α : A ⟶ A) = P.Base (𝟙 A) from hα.1,
      P.Base_id, Category.comp_id, show P.Base ζn = P.Base (𝟙 A) from hζb,
      P.Base_id, Category.comp_id] at h
    rw [P.Base_id]
    simpa using h.symm
  have hβm : β ∈ OTri P A := ⟨hβb, hβs.1⟩
  -- Frobenius-normalized で `ζ n ≫ β^n = β ≫ ζ n`
  have hfnb := hfn ζn hζb β hβm
  rw [hζd] at hfnb
  refine ⟨β, ⟨hβm, ?_⟩, ?_⟩
  · -- 全射性: `ζ n` が epi
    refine (cancel_epi (ζn : A ⟶ A)).mp ?_
    rw [hfnb, hsq]
  · -- 単射性: (b) の `∃!`
    rintro γ ⟨hγm, hγp⟩
    refine huniq γ ⟨⟨hγm.2, ?_⟩, ?_⟩
    · show IsIso (P.Base (γ : A ⟶ A))
      rw [show P.Base (γ : A ⟶ A) = P.Base (𝟙 A) from hγm.1, P.Base_id]
      infer_instance
    · have hfng := hfn ζn hζb γ hγm
      rw [hζd, hγp] at hfng
      exact hfng

include P in
/-- ★★★**`Proposition 1.10, (iii)` の moreover 完成** —— `𝒪^×(A)` も perfect。

★`𝒪^▷(A)` の結果に、★**「`β^n` が同型なら `β` も同型」**(`isIso_of_pow_isIso`)を
足すだけ。★**単元の側は `𝒪^▷` の側から自動で従う。** -/
theorem prop_1_10_iii_otimes_perfect (F : FrobenioidCore P) {A : C}
    (hiso : ∀ X : C, IsIsotropic P X)
    (hperf : IsPerfectObj P A) (hfn : IsFrobeniusNormalized P A)
    (n : ℕ+) (ζn : End A) (hζb : IsBaseIdentity P ζn) (hζf : IsFrobeniusType P ζn)
    (hζd : P.degFr ζn = n) (α : End A) (hα : α ∈ OTimes P A) :
    ∃! β : End A, β ∈ OTimes P A ∧ (β ^ (n : ℕ) : End A) = α := by
  obtain ⟨β, ⟨hβm, hβp⟩, huniq⟩ :=
    prop_1_10_iii_otri_perfect P hperf hfn n ζn hζb hζf hζd α hα.1
  -- `β^n = α` は単元なので `β` も単元
  haveI hpow : IsIso ((β ^ (n : ℕ) : End A) : A ⟶ A) := by
    rw [hβp]
    exact (CategoryTheory.isUnit_iff_isIso (α : End A)).mp hα.2
  have hβiso : IsIso ((β : End A) : A ⟶ A) :=
    isIso_of_pow_isIso P F hiso β hβm.1 hβm.2 n.pos hpow
  refine ⟨β, ⟨⟨hβm, (CategoryTheory.isUnit_iff_isIso (β : End A)).mpr hβiso⟩, hβp⟩, ?_⟩
  rintro γ ⟨hγm, hγp⟩
  exact huniq γ ⟨hγm.1, hγp⟩

/-! ### ★★(iii) の還元 —— 原文「we may assume that A is Frobenius-trivial」

原文 (FrdI p.35):
> of Frobenius-trivial objects [cf. Definition 1.3, (i), (a), (b); the isomorphism of Def-

原文 (FrdI p.35):
> inition 1.3, (iii), (c)], we may assume that A is Frobenius-trivial. Now the fact that

★★**取り下げの原因の1つを埋める**(2026-08-16)。
`prop_1_10_iii_otri_perfect` は「次数 `n` の base-identity Frobenius 型自己射 `ζ n` が
ある」を**仮定に置いていた** —— それは `A` が Frobenius-trivial だということである。
★原文はそこを**還元して**いる。その還元をここで実装する。

★**原文が名指す 3 つが、そのまま 3 段になる**:
1. `Definition 1.3, (i), (a)`(`baseSurj`) —— `Base(A)` を底に持つ
   **Frobenius-trivial な `A₀`** がある
2. `Definition 1.3, (i), (b)`(`preStepSpan`) —— その底の同型は
   **pre-step の span** `A₀ ← X → A` に持ち上がる
3. `Definition 1.3, (iii), (c)`(`otriFwd` / `otriBwd`) —— **co-angular** pre-step は
   `𝒪^▷` の全単射を与える

★★**`isotropic` の仮定がここで効く**。原文の角括弧
「so all morphisms of C are co-angular — cf. Proposition 1.4, (i)」がそれである ——
span の pre-step が co-angular だと分かって初めて 3 が使える。
★**原文は「isotropic」を moreover の仮定として置いているが、
その最初の用途はこの還元である**(perfect 性の議論そのものではない)。

★★**「一意存在が移る」ことだけが要る**。全単射そのもの(`MulEquiv`)は要らない ——
`otriFwd` / `otriBwd` の 2 つの `∃!` と、四角形が積と可換すること
(`φ ≫ (β₂ ≫ β₁) = (α₂ ≫ α₁) ≫ φ`)があれば、`n` 乗根の一意存在は渡る。
-/

section PowTransfer

variable {M N : Type*} [Monoid M] [Monoid N]

/-- ★**全単射の graph に沿って「`n` 乗根の一意存在」が移る**。

★`R` は `Definition 1.3, (iii), (c)` の四角形「`φ ≫ β = α ≫ φ`」であり、
`otriFwd` / `otriBwd` がちょうど**両向きの一意性**を与える。
★積と可換する(`hmul`)ことから `n` 乗とも可換し、両側の一意性で挟める。

★**`MulEquiv` を作らないのは、`R` が `∃!` の形でしか与えられないからである** ——
写像を取り出すには選択が要るが、主張は選択なしで書ける。 -/
theorem exists_unique_pow_transfer (R : M → N → Prop)
    (hfwd : ∀ a : M, ∃! b : N, R a b) (hbwd : ∀ b : N, ∃! a : M, R a b)
    (hone : R 1 1)
    (hmul : ∀ (a₁ a₂ : M) (b₁ b₂ : N), R a₁ b₁ → R a₂ b₂ → R (a₁ * a₂) (b₁ * b₂))
    (n : ℕ) (h : ∀ a : M, ∃! x : M, x ^ n = a) (b : N) : ∃! y : N, y ^ n = b := by
  have hpow : ∀ (a : M) (b : N), R a b → ∀ k : ℕ, R (a ^ k) (b ^ k) := by
    intro a b hab k
    induction k with
    | zero => rw [pow_zero, pow_zero]; exact hone
    | succ k ih => rw [pow_succ, pow_succ]; exact hmul _ _ _ _ ih hab
  obtain ⟨a, hab, hau⟩ := hbwd b
  obtain ⟨x, hxa, hxu⟩ := h a
  obtain ⟨y, hxy, hyu⟩ := hfwd x
  obtain ⟨z, -, hzu⟩ := hfwd (x ^ n)
  refine ⟨y, ?_, ?_⟩
  · show y ^ n = b
    rw [hzu _ (hpow x y hxy n), hzu b (by rw [hxa]; exact hab)]
  · intro y' hy'
    have hy'' : y' ^ n = b := hy'
    obtain ⟨a', ha'y, -⟩ := hbwd y'
    have ha'n : a' ^ n = a := hau _ (by rw [← hy'']; exact hpow a' y' ha'y n)
    exact hyu y' (by rw [← hxu a' ha'n]; exact ha'y)

end PowTransfer

/-- ★`𝒪^▷` の「`n` 乗根の一意存在」の 2 つの綴りが同値であること ——
`End A` の元＋所属で書いた形と、部分モノイドの元で書いた形。

★**前者は原文の言い回しに近く、後者は移送の補題に渡しやすい。** -/
theorem otri_existsUnique_pow_iff {A : C} (n : ℕ) :
    (∀ α ∈ OTri P A, ∃! β : End A, β ∈ OTri P A ∧ (β ^ n : End A) = α)
      ↔ (∀ a : OTri P A, ∃! x : OTri P A, x ^ n = a) := by
  constructor
  · intro h a
    obtain ⟨β, ⟨hβm, hβp⟩, hβu⟩ := h (a : End A) a.2
    refine ⟨⟨β, hβm⟩, Subtype.ext (by simpa using hβp), fun y hy => ?_⟩
    exact Subtype.ext (hβu (y : End A) ⟨y.2, by simpa using congrArg Subtype.val hy⟩)
  · intro h α hα
    obtain ⟨x, hx, hxu⟩ := h ⟨α, hα⟩
    refine ⟨(x : End A), ⟨x.2, by simpa using congrArg Subtype.val hx⟩, ?_⟩
    rintro β ⟨hβm, hβp⟩
    exact congrArg Subtype.val (hxu ⟨β, hβm⟩ (Subtype.ext (by simpa using hβp)))

/-- ★★**`Definition 1.3, (iii), (c)` を「perfect 性の移送」として使う**。

co-angular pre-step `φ : A ⟶ B` に沿って、`𝒪^▷` の `n` 乗根の一意存在は
★**両向きに**移る(`otriFwd` と `otriBwd` の両方があるため)。 -/
theorem otri_pow_transfer (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hc : IsCoAngular P φ) (hs : IsPreStep P φ) (n : ℕ) :
    (∀ a : OTri P A, ∃! x : OTri P A, x ^ n = a) ↔
      (∀ b : OTri P B, ∃! y : OTri P B, y ^ n = b) := by
  set R : OTri P A → OTri P B → Prop :=
    fun a b => (φ ≫ (b : End B) : A ⟶ B) = (a : End A) ≫ φ with hR
  have hfwd : ∀ a : OTri P A, ∃! b : OTri P B, R a b := by
    intro a
    obtain ⟨b, ⟨hbm, hbr⟩, hbu⟩ := F.otriFwd φ hc hs (a : End A) a.2
    exact ⟨⟨b, hbm⟩, hbr, fun y hy => Subtype.ext (hbu (y : End B) ⟨y.2, hy⟩)⟩
  have hbwd : ∀ b : OTri P B, ∃! a : OTri P A, R a b := by
    intro b
    obtain ⟨a, ⟨ham, har⟩, hau⟩ := F.otriBwd φ hc hs (b : End B) b.2
    exact ⟨⟨a, ham⟩, har, fun y hy => Subtype.ext (hau (y : End A) ⟨y.2, hy⟩)⟩
  have hone : R 1 1 := by
    show (φ ≫ (𝟙 B : End B) : A ⟶ B) = (𝟙 A : End A) ≫ φ
    simp
  have hmul : ∀ (a₁ a₂ : OTri P A) (b₁ b₂ : OTri P B),
      R a₁ b₁ → R a₂ b₂ → R (a₁ * a₂) (b₁ * b₂) := by
    intro a₁ a₂ b₁ b₂ h₁ h₂
    show (φ ≫ ((b₂ : End B) ≫ (b₁ : End B)) : A ⟶ B)
      = ((a₂ : End A) ≫ (a₁ : End A)) ≫ φ
    rw [← Category.assoc φ, show (φ ≫ (b₂ : End B) : A ⟶ B) = (a₂ : End A) ≫ φ from h₂,
      Category.assoc, show (φ ≫ (b₁ : End B) : A ⟶ B) = (a₁ : End A) ≫ φ from h₁,
      ← Category.assoc]
  exact ⟨fun h => exists_unique_pow_transfer R hfwd hbwd hone hmul n h,
    fun h => exists_unique_pow_transfer (fun b a => R a b) hbwd hfwd hone
      (fun b₁ b₂ a₁ a₂ h₁ h₂ => hmul a₁ a₂ b₁ b₂ h₁ h₂) n h⟩

/-- ★★★**`Proposition 1.10, (iii)` の moreover —— 還元まで込めた完成形**(`𝒪^▷`)。

★原文の仮定そのもの(「`𝒞` が perfect・isotropic・Frobenius-normalized 型」)から、
★**任意の対象 `A`** について `𝒪^▷(A)` で `n` 乗根が一意に存在する。
`ζ n` の存在はもはや仮定ではなく、`Definition 1.3, (i), (a)` から**導かれる**。 -/
theorem prop_1_10_iii_otri_perfect_of_type (F : FrobenioidCore P)
    (hpt : IsOfPerfectType P) (hiso : ∀ X : C, IsIsotropic P X)
    (hfnt : ∀ X : C, IsFrobeniusNormalized P X) (A : C) (n : ℕ+)
    (α : End A) (hα : α ∈ OTri P A) :
    ∃! β : End A, β ∈ OTri P A ∧ (β ^ (n : ℕ) : End A) = α := by
  -- 1 段: `Definition 1.3, (i), (a)` —— 同じ底を持つ Frobenius-trivial な `A₀`
  obtain ⟨A₀, ⟨ζ, hζd, hζbf⟩, ⟨e⟩⟩ := F.baseSurj (P.toElem.obj A).base
  have hA₀ : ∀ a : OTri P A₀, ∃! x : OTri P A₀, x ^ (n : ℕ) = a :=
    (otri_existsUnique_pow_iff P (n : ℕ)).mp (fun β hβ =>
      prop_1_10_iii_otri_perfect P (hpt A₀) (hfnt A₀) n (ζ n)
        (hζbf n).1 (hζbf n).2 (hζd n) β hβ)
  -- 2 段: `Definition 1.3, (i), (b)` —— pre-step の span `A₀ ← X → A`
  obtain ⟨X, φ, ψ, hφs, hψs, -⟩ := F.preStepSpan A₀ A e.hom inferInstance
  -- 3 段: isotropic ⟹ co-angular(`Proposition 1.4, (i)`)、`(iii)(c)` で両向きに移送
  have hφc : IsCoAngular P φ := prop_1_4_i P φ (fun Y _ => hiso Y)
  have hψc : IsCoAngular P ψ := prop_1_4_i P ψ (fun Y _ => hiso Y)
  exact (otri_existsUnique_pow_iff P (n : ℕ)).mpr
    ((otri_pow_transfer P F ψ hψc hψs (n : ℕ)).mp
      ((otri_pow_transfer P F φ hφc hφs (n : ℕ)).mpr hA₀)) α hα

/-- ★★★**`Proposition 1.10, (iii)` の moreover —— 還元まで込めた完成形**(`𝒪^×`)。

★`𝒪^▷` の側に「`β^n` が同型なら `β` も同型」を足すだけ。 -/
theorem prop_1_10_iii_otimes_perfect_of_type (F : FrobenioidCore P)
    (hpt : IsOfPerfectType P) (hiso : ∀ X : C, IsIsotropic P X)
    (hfnt : ∀ X : C, IsFrobeniusNormalized P X) (A : C) (n : ℕ+)
    (α : End A) (hα : α ∈ OTimes P A) :
    ∃! β : End A, β ∈ OTimes P A ∧ (β ^ (n : ℕ) : End A) = α := by
  obtain ⟨β, ⟨hβm, hβp⟩, huniq⟩ :=
    prop_1_10_iii_otri_perfect_of_type P F hpt hiso hfnt A n α hα.1
  haveI hpow : IsIso ((β ^ (n : ℕ) : End A) : A ⟶ A) := by
    rw [hβp]
    exact (CategoryTheory.isUnit_iff_isIso (α : End A)).mp hα.2
  have hβiso : IsIso ((β : End A) : A ⟶ A) :=
    isIso_of_pow_isIso P F hiso β hβm.1 hβm.2 n.pos hpow
  refine ⟨β, ⟨⟨hβm, (CategoryTheory.isUnit_iff_isIso (β : End A)).mpr hβiso⟩, hβp⟩, ?_⟩
  rintro γ ⟨hγm, hγp⟩
  exact huniq γ ⟨hγm.1, hγp⟩

end ABC3.Found.FrdI
