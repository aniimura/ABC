/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.RatTower
import Mathlib.Algebra.Algebra.Subalgebra.Directed

/-!
# [GenEll] Remark 1.5.1 —— **余極限へ向けた整合性**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

## ★★★ここまでで建ったもの

`RatTower.lean` の図式 `ratTowerDiagram`（`n ↦ ℤ[1/n!]`）と余錐 `ratTowerCocone`
（頂点 `ℚ`）に対し、**任意の余錐からの写像**を作る。

| 主張 | 宣言 |
|---|---|
| 図式の射は「`ℚ` の中では恒等」 | ★`ratTowerDiagram_map_coe` |
| 段を上げても値は変わらない | ★`cocone_app_le` |
| 同じ有理数を表す元は同じ値へ | ★★`cocone_app_compat` |
| ★**余極限からの写像** | ★★★`descFun` |
| どの段で読んでも同じ | ★★`descFun_eq` |
| `ℤ`-部分代数版（`iSupLift` へ渡す形） | `ratTowerAlg` ほか |

★★`descFun` が**環準同型であること**が残っている。★機構は分かっている——
`exists_ratTower_of_finset` で `a`・`b`・`a+b` を共通の段に集め、
その段の射が環準同型であることを使う。

## ★★★★★配管の記録（`tools/lean-idioms.md` にも 1 行足した）

★**`ratTower n = Subring.closure {...}` なので `CommRingCat.of (ratTower n)` の
インスタンスが核にとって重い。** 2026-08-27 に 2 回落ちた:

1. `Functor` の構造体を素朴に書く（`map_id` / `map_comp` を `ext; rfl`）
   —— 核判定が 100 秒を超えて落ちない。★`Functor.ofSequence` に替えたら 0.07 秒。
2. `towerMk n q hq : ratTowerDiagram.obj n := ⟨q, hq⟩` という補助を挟む
   —— こちらも核判定が 95 秒。★`show ratTower n from ⟨q, hq⟩` と直に書けば 0.05 秒。

★★★**教訓: `Subring.closure` を台にした型を圏論の構造に載せるときは、
補助定義を挟まず、既にある構成子（`ofSequence` 等）に寄せる。**
★根本的に直すなら `ratTower` を明示的な `carrier` で定義し直すことになる。
-/

namespace ABC3.Found.GenEll

open CategoryTheory Limits

/-! ## ★★`ℤ`-部分代数版 —— `Subalgebra.iSupLift` へ渡すため -/

/-- ★`ℤ[1/n!]` を `ℤ`-部分代数として。

★`Subalgebra.iSupLift`（`Mathlib/Algebra/Algebra/Subalgebra/Directed.lean`）は
**`Subalgebra R A` 用で `Subring` 版が無い**ので、橋を渡しておく。 -/
def ratTowerAlg (n : ℕ) : Subalgebra ℤ ℚ := subalgebraOfSubring (ratTower n)

theorem mem_ratTowerAlg {n : ℕ} {q : ℚ} : q ∈ ratTowerAlg n ↔ q ∈ ratTower n :=
  mem_subalgebraOfSubring

theorem ratTowerAlg_mono {n m : ℕ} (h : n ≤ m) : ratTowerAlg n ≤ ratTowerAlg m := by
  intro x hx
  exact mem_ratTowerAlg.2 (ratTower_mono h (mem_ratTowerAlg.1 hx))

theorem ratTowerAlg_directed : Directed (· ≤ ·) ratTowerAlg := by
  intro n m
  exact ⟨max n m, ratTowerAlg_mono (le_max_left n m), ratTowerAlg_mono (le_max_right n m)⟩

theorem iSup_ratTowerAlg_eq_top : (⨆ n : ℕ, ratTowerAlg n) = ⊤ := by
  refine eq_top_iff.2 (fun q _ => ?_)
  exact le_iSup ratTowerAlg q.den (mem_ratTowerAlg.2 (mem_ratTower_den q))

/-! ## ★★★★図式の射は「`ℚ` の中では恒等」 -/

/-- ★★**図式の射は `ℚ` の中では恒等**——余錐の自然性から読む。

★`Functor.ofSequence` で作った射は「連続する包含の合成」なので、
直に `Subring.inclusion` と比べるには帰納法が要る。★★**余錐の自然性を経由すれば要らない。** -/
theorem ratTowerDiagram_map_coe {n m : ℕ} (h : n ⟶ m) (x : ratTower n) :
    ((show (ratTower m) from (ratTowerDiagram.map h).hom x) : ℚ) = (x : ℚ) := by
  have hnat := ratTowerCocone.ι.naturality h
  have h2 := congrArg (fun (f : CommRingCat.of (ratTower n) ⟶ CommRingCat.of ℚ) => f.hom x) hnat
  simp only [ratTowerCocone, Function.comp_apply, Category.comp_id,
    Functor.const_obj_map] at h2
  exact h2

/-! ## ★★★★★★任意の余錐との整合性 -/

/-- ★★**段を上げても値は変わらない**。 -/
theorem cocone_app_le {s : Cocone ratTowerDiagram} {i j : ℕ} (hij : i ≤ j)
    (x : ratTowerDiagram.obj i) :
    (s.ι.app i).hom x = (s.ι.app j).hom ((ratTowerDiagram.map (homOfLE hij)).hom x) :=
  (congrArg (fun f => CommRingCat.Hom.hom f x) (s.w (homOfLE hij))).symm

/-- ★★★**同じ有理数を表す 2 つの段の元は同じ値に送られる**。

★★これが「余極限からの写像が well-defined」の中身である。
★機構は有向性——`max i j` へ両方を上げて `ratTowerDiagram_map_coe` で比べる。 -/
theorem cocone_app_compat {s : Cocone ratTowerDiagram} {i j : ℕ}
    (x : ratTowerDiagram.obj i) (y : ratTowerDiagram.obj j)
    (hxy : ((show ratTower i from x) : ℚ) = ((show ratTower j from y) : ℚ)) :
    (s.ι.app i).hom x = (s.ι.app j).hom y := by
  have hik : i ≤ max i j := le_max_left i j
  have hjk : j ≤ max i j := le_max_right i j
  have hx := cocone_app_le (s := s) hik x
  have hy := cocone_app_le (s := s) hjk y
  rw [hx, hy]
  congr 1
  apply Subtype.ext
  rw [ratTowerDiagram_map_coe, ratTowerDiagram_map_coe]
  exact hxy

/-! ## ★★★★★★★★余極限からの写像 -/

/-- ★★★★**余極限からの写像** —— `q` をその分母の段で読む。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX (Q) depends only on the pair

★`mem_ratTower_den`（どの有理数も `ℤ[1/q.den!]` に入る）で段が決まる。 -/
noncomputable def descFun (s : Cocone ratTowerDiagram) (q : ℚ) : s.pt :=
  (s.ι.app q.den).hom (show ratTower q.den from ⟨q, mem_ratTower_den q⟩)

/-- ★★**どの段で読んでも同じ** —— `cocone_app_compat` そのもの。 -/
theorem descFun_eq (s : Cocone ratTowerDiagram) (n : ℕ) (q : ℚ) (hq : q ∈ ratTower n) :
    descFun s q = (s.ι.app n).hom (show ratTower n from ⟨q, hq⟩) :=
  cocone_app_compat (s := s) (show ratTowerDiagram.obj q.den from ⟨q, mem_ratTower_den q⟩)
    (show ratTowerDiagram.obj n from ⟨q, hq⟩) rfl

/-! ### ★出典の紐付け(`.src`)

★★**項目全体の `.src` は置かない。** 残っているのは
(a) `descFun` が環準同型であること、(b) `IsColimit` に束ねること、
(c) `Spec` で `IsLimit` に移すこと、(d) `X ×_ℤ (−)` で `X_ℚ` の図式にすること、
(e) mathlib の `Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation` を当てること、
(f) 因子 `D` の降下、(g) `Σ` の外での conductor の一致、である。 -/

def ratTowerAlg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(ℤ[1/n!] の ℤ-部分代数版——Subalgebra.iSupLift へ渡すため)",
    sectionId := "genell-rem-1-5-1" }

def ratTowerDiagram_map_coe.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(図式の射は ℚ の中では恒等)",
    sectionId := "genell-rem-1-5-1" }

def cocone_app_compat.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(同じ有理数を表す元は同じ値へ——well-defined 性)",
    sectionId := "genell-rem-1-5-1" }

def descFun.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(余極限からの写像——環準同型であることは含まない)",
    sectionId := "genell-rem-1-5-1" }

def descFun.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "cocone_app_compat(同じ有理数を表す元は同じ値へ)"
      (.inProject "ABC3" "ABC3.Found.GenEll.cocone_app_compat") 9,
    .citation "[ABC3]" "mem_ratTower_den(どの有理数も ℤ[1/q.den!] に入る)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mem_ratTower_den") 9,
    .implicitStep
      ("★★descFun が環準同型であることが残っている。機構は分かっている——" ++
       "exists_ratTower_of_finset で a・b・a+b を共通の段に集め、" ++
       "その段の射が環準同型であることを使う") 9,
    .implicitStep
      ("★★★★★配管の記録: ratTower n = Subring.closure {...} なので " ++
       "CommRingCat.of (ratTower n) のインスタンスが核にとって重い。" ++
       "2026-08-27 に 2 回落ちた(Functor の構造体を素朴に書く / towerMk の補助を挟む)。" ++
       "既にある構成子(Functor.ofSequence 等)に寄せ、補助定義を挟まないこと") 9,
    .implicitStep
      ("★残りの段: (a) 環準同型性 (b) IsColimit (c) Spec で IsLimit へ " ++
       "(d) X ×_ℤ (−) で X_ℚ の図式へ (e) mathlib の " ++
       "Scheme.exists_π_app_comp_eq_of_locallyOfFinitePresentation " ++
       "(f) 因子 D の降下 (g) Σ の外での conductor の一致") 9 ]

end ABC3.Found.GenEll
