import ABC3.Interface.GaloisRep.Reduction

/-!
# 退化封じの検査 —— **G7 の `ω_E` は曲線と結ばれていない**(`Check`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★2026-08-18 に見つかった 3 つ目の穴(同じ型)

`SemistableModelData` の `omega`(Néron 微分)を縛っているのは
**`omega_rank_one`(階数 1)だけ**である。

★★★**`omega L W := 𝓞 L`(曲線を無視する定数)がこれを満たす。**
★したがって `ω_E` は**曲線と一切結ばれていない**。

## ★★これで 3 回目である

| # | 場所 | 穴 | 塞いだ日 |
|---|---|---|---|
| 1 | C1 `evalAffine` | 値がアフィンでしか固定されていなかった | 2026-08-17 |
| 2 | B1 `Pic` | 非アフィンで自由だった(`sheafOf` を追加) | 2026-08-17 |
| 3 | B1 `pullback` | `sheafOf` と結ばれていなかった | 2026-08-18 |
| 4 | **G7 `omega`** | **曲線と結ばれていない**(本ファイル) | ★★2026-08-26 |

## ★★★★★★一般規則(この 4 例から抽出)

> **`→ Type` の posit(台を仮定する欄)は、必ず
> 「その欄と入力データの両方を言及する条件」を 1 本以上持たねばならない。**

★階数・濃度・群構造だけでは足りない——それらは**入力を無視する定数**で満たせる。

## ★注意(限界の明示)

★本ファイルは「`omega` の欄が入力と結ばれていない」ことを示すのであって、
「`SemistableModelData` 全体が退化 witness を持つ」ことではない
(それには `TateCurveData` の witness が要るが、G6 は未実装)。
★★`Check/Arakelov/PullbackNondegenerate.lean` と同じ立場である。

## ★★★★★★2026-08-26:穴は塞いだ

`SemistableModelData` の欄を **`omegaFrac : … → Submodule (𝓞 L) L`**(分数イデアル)に
差し替え、**変数変換での変わり方**

    x ∈ ω_{C•E}  ↔  u·x ∈ ω_E

を課した。★★★★これは Néron 微分について**真**である(不変微分が `ω' = u·ω` と変わるから)。
★★本ファイル後半で、**この 3 条件を満たす定数 witness は存在しない**ことを証明する
(`no_constant_omegaFrac`)——規則「`→ Type` の posit は入力と結ぶ条件を持て」の 4 例目が閉じた。
-/

namespace ABC3.Check.GaloisRep

open ABC3.Meta WeierstrassCurve NumberField

/-- ★`SemistableModelData` から `omega` の 4 欄だけを取り出した構造。 -/
structure WeakOmega where
  /-- Néron 微分の台。 -/
  omega : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → Type
  /-- 加法群の構造。 -/
  omegaAddCommGroup : (L : Type) → [Field L] → [NumberField L] →
    (W : WeierstrassCurve L) → AddCommGroup (omega L W)
  /-- `𝓞_L` 加群の構造。 -/
  omegaModule : (L : Type) → [Field L] → [NumberField L] →
    (W : WeierstrassCurve L) →
    @Module (𝓞 L) (omega L W) _ (omegaAddCommGroup L W).toAddCommMonoid
  /-- ★階数 1。 -/
  omega_rank_one : ∀ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L),
    @Module.rank (𝓞 L) (omega L W) _ (omegaAddCommGroup L W).toAddCommMonoid
      (omegaModule L W) = 1

/-- ★★★★**曲線を無視する witness** —— `ω_E := 𝓞_L`。

★★★階数 1 の条件は満たすが、**曲線 `W` を一切見ていない**。 -/
def constOmega : WeakOmega where
  omega L _ _ _ := 𝓞 L
  omegaAddCommGroup _ _ _ _ := inferInstance
  omegaModule _ _ _ _ := inferInstance
  omega_rank_one _ _ _ _ := Module.rank_self _

/-- ★★★★★**穴は実在した** —— `omega` は曲線に依存しない witness を通す。

★これが「`→ Type` の posit は入力と結ぶ条件を 1 本以上持て」という規則の 4 例目である。 -/
theorem omega_nondegenerate_hole :
    ∃ W : WeakOmega, ∀ (L : Type) [Field L] [NumberField L]
      (E E' : WeierstrassCurve L), W.omega L E = W.omega L E' :=
  ⟨constOmega, by intros; rfl⟩

/-! ## ★★★★★★2026-08-26:塞いだ形の検査 -/

section Strong

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★**変数変換則は「全体」を強制する**——定数 witness の逃げ場を奪う鍵。 -/
theorem law_forces_top (I : Submodule (𝓞 L) L) (hne : I ≠ ⊥)
    (hlaw : ∀ (u : Lˣ) (x : L), x ∈ I ↔ (u : L) * x ∈ I) : I = ⊤ := by
  obtain ⟨c, hcI, hc0⟩ := Submodule.exists_mem_ne_zero_of_ne_bot hne
  refine Submodule.eq_top_iff'.2 fun y => ?_
  rcases eq_or_ne y 0 with rfl | hy
  · exact I.zero_mem
  · have hu : IsUnit (c / y) := (isUnit_iff_ne_zero).2 (div_ne_zero hc0 hy)
    refine (hlaw hu.unit y).2 ?_
    rw [hu.unit_spec, div_mul_cancel₀ _ hy]
    exact hcI

/-- ★★★★`L` は `𝓞_L` 上有限生成でない——`⊤` は分数イデアルになれない。 -/
theorem not_fg_top : ¬ (⊤ : Submodule (𝓞 L) L).FG := by
  intro hfg
  haveI : Module.Finite (𝓞 L) L := ⟨hfg⟩
  have hsurj : Function.Surjective (algebraMap (𝓞 L) L) := by
    intro y
    have hint : IsIntegral (𝓞 L) y := IsIntegral.of_finite _ _
    exact IsIntegrallyClosed.isIntegral_iff.1 hint
  have hbij : Function.Bijective (algebraMap (𝓞 L) L) :=
    ⟨IsFractionRing.injective _ _, hsurj⟩
  exact NumberField.RingOfIntegers.not_isField L
    (MulEquiv.isField (Field.toIsField L) (RingEquiv.ofBijective _ hbij).toMulEquiv)

/-- ★★★★★★**定数 witness はもう通らない**——穴は塞がった。

★曲線を無視する分数イデアル `I` は、変数変換則を満たすと `⊤` になり、
有限生成でなくなる。 -/
theorem no_constant_omegaFrac (I : Submodule (𝓞 L) L) (hfg : I.FG) (hne : I ≠ ⊥) :
    ¬ (∀ (W : WeierstrassCurve L) (C : WeierstrassCurve.VariableChange L) (x : L),
        x ∈ I ↔ ((C.u : L) * x) ∈ I) := by
  intro hlaw
  have hlaw' : ∀ (u : Lˣ) (x : L), x ∈ I ↔ (u : L) * x ∈ I := by
    intro u x
    exact hlaw ⟨0, 0, 0, 0, 0⟩ ⟨u, 0, 0, 0⟩ x
  exact not_fg_top (law_forces_top I hne hlaw' ▸ hfg)

end Strong

/-! ## ★出典の紐付け(`.src`) -/

def omega_nondegenerate_hole.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(層 G——Néron 微分が曲線と結ばれていないこと)",
    sectionId := "genell-def-3-3" }

def no_constant_omegaFrac.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Néron 微分の欄を曲線に結びつけたこと——穴を塞いだ検査)",
    sectionId := "genell-def-3-3" }

end ABC3.Check.GaloisRep
